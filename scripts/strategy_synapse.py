import time
from tqdm import tqdm
import argparse
import numpy as np
import pandas as pd
from caveclient import CAVEclient

parser = argparse.ArgumentParser()
parser.add_argument("--ver", "-v", help="materialization version",type=int, default=1412)
parser.add_argument("--unpf_samp", "-n", help="number of unproofread cells to sample",type=int, default=1000)
parser.add_argument("--seed", "-s", help="RNG seed",type=int, default=0)
args = vars(parser.parse_args())
ver = int(args["ver"])
unpf_samp = int(args["unpf_samp"])
seed = int(args["seed"])

client = CAVEclient("minnie65_public")
client.materialize.get_versions()

client.version = ver

# load dataframe of all cells with proofread axons/dendrites
axon_pf_df = client.materialize.tables.proofreading_status_and_strategy(status_axon="t").query() \
    [["pt_root_id","strategy_axon"]].rename(columns={"strategy_axon":"strategy"})
ddrt_pf_df = client.materialize.tables.proofreading_status_and_strategy(status_dendrite="t").query() \
    [["pt_root_id","strategy_dendrite"]].rename(columns={"strategy_dendrite":"strategy"})

# load dataframe of cell types
type_df = client.materialize.tables.aibs_metamodel_mtypes_v661_v2().query() \
    [["pt_root_id","cell_type"]]
type_df.drop_duplicates(subset="pt_root_id", inplace=True)
non_neur_type_df = client.materialize.tables.aibs_metamodel_celltypes_v661(classification_system="nonneuron") \
    .query()[["pt_root_id","cell_type"]]
non_neur_type_df["cell_type"] = "nonneuron"
type_df = pd.concat([type_df, non_neur_type_df], ignore_index=True)
axon_pf_df = axon_pf_df.merge(type_df, on="pt_root_id", how="left")
ddrt_pf_df = ddrt_pf_df.merge(type_df, on="pt_root_id", how="left")

# load and sample dataframe of unproofread neurons
unpf_df = client.materialize.tables.nucleus_detection_v0().query()[["pt_root_id"]]
unpf_df = unpf_df[unpf_df["pt_root_id"] != 0]
unpf_df = unpf_df.merge(type_df, on="pt_root_id", how="left")
unpf_df = unpf_df[unpf_df["cell_type"] != "nonneuron"]
unpf_df = unpf_df[~unpf_df["pt_root_id"].isin(axon_pf_df["pt_root_id"]) &
                  ~unpf_df["pt_root_id"].isin(ddrt_pf_df["pt_root_id"])]
unpf_df = unpf_df.iloc[np.random.default_rng(seed).choice(len(unpf_df), size=unpf_samp, replace=False)]
unpf_df["strategy"] = "unproofread"

# concatenate unproofread cells to axon and dendrite proofread dataframes
axon_pf_df = pd.concat([axon_pf_df, unpf_df], ignore_index=True)
ddrt_pf_df = pd.concat([ddrt_pf_df, unpf_df], ignore_index=True)

# fill NaN values in cell_type columns
axon_pf_df["cell_type"] = axon_pf_df["cell_type"].fillna("Unknown")
ddrt_pf_df["cell_type"] = ddrt_pf_df["cell_type"].fillna("Unknown")

# assign Exc/Inh labels based on cell type
for df in [axon_pf_df, ddrt_pf_df]:
    df["Exc/Inh"] = ""
    exc_cells = df["cell_type"].str.contains("L").values.astype(bool)
    df.loc[exc_cells,"Exc/Inh"] = "Exc"
    inh_cells = df["cell_type"].str.contains("TC").values.astype(bool)
    df.loc[inh_cells,"Exc/Inh"] = "Inh"
    unk_cells = df["cell_type"].values == "Unknown"
    df.loc[unk_cells,"Exc/Inh"] = "Unknown"
    
# compute number of input/output synapses for each cell
axon_pf_syns = np.zeros_like(axon_pf_df["pt_root_id"].values)
ddrt_pf_syns = np.zeros_like(ddrt_pf_df["pt_root_id"].values)

for i,pt_root_id in tqdm(enumerate(axon_pf_df["pt_root_id"].values)):
    axon_pf_syns[i] = len(client.materialize.synapse_query(pre_ids=pt_root_id, remove_autapses=True))
for i,pt_root_id in tqdm(enumerate(ddrt_pf_df["pt_root_id"].values)):
    ddrt_pf_syns[i] = len(client.materialize.synapse_query(post_ids=pt_root_id, remove_autapses=True))
    
axon_pf_df["num_synapses"] = axon_pf_syns
ddrt_pf_df["num_synapses"] = ddrt_pf_syns

axon_pf_df.to_csv("./results/pf_strat_axon.csv", index=False)
ddrt_pf_df.to_csv("./results/pf_strat_ddrt.csv", index=False)
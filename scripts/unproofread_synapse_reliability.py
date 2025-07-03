import pickle
import time
from tqdm import tqdm
import argparse
import numpy as np
import pandas as pd
from caveclient import CAVEclient

parser = argparse.ArgumentParser()
parser.add_argument("--curr_ver", "-c", help="current materialization version to compare to",type=int, default=1412)
parser.add_argument("--prev_ver", "-p", help="previous materialization version to compare to",type=int, default=1412)
parser.add_argument("--num_samp", "-n", help="number of presynaptic cells to sample",type=int, default=100)
parser.add_argument("--num_updt", "-nu", help="number of cell ids to update",type=int, default=2000)
parser.add_argument("--seed", "-s", help="RNG seed",type=int, default=0)
args = vars(parser.parse_args())
curr_ver = int(args["curr_ver"])
prev_ver = int(args["prev_ver"])
num_samp = int(args["num_samp"])
num_updt = int(args["num_updt"])
seed = int(args["seed"])

client = CAVEclient("minnie65_public")
client.materialize.get_versions()

client.version = curr_ver

# load dataframe of cell ids in cortical column
column_df = client.materialize.tables.allen_column_mtypes_v2().query()[["pt_root_id"]]

# load dataframe of cell types
type_df = client.materialize.tables.aibs_metamodel_mtypes_v661_v2() \
    .query(select_columns={"nucleus_detection_v0": ["pt_root_id"],
                           "aibs_metamodel_mtypes_v661_v2": ["cell_type"]})
type_df.drop_duplicates(subset="pt_root_id", inplace=True)

# load previous and current dataframes of cell ids with unproofread axons
curr_proof_df = client.materialize.tables.proofreading_status_and_strategy(status_axon="t") \
    .query()[["pt_root_id"]]
curr_unpf_df = client.materialize.tables.nucleus_detection_v0() \
    .query()[["pt_root_id"]]
curr_unpf_df = curr_unpf_df[~curr_unpf_df["pt_root_id"].isin(curr_proof_df["pt_root_id"])]
if prev_ver >= 1078:
    prev_proof_df = client.materialize.tables.proofreading_status_and_strategy(status_axon="t") \
        .query(timestamp=client.materialize.get_version_metadata(prev_ver)["time_stamp"]) \
            [["pt_root_id"]]
else:
    client.version = prev_ver
    prev_proof_df = client.materialize.query_table('proofreading_status_public_release',
                                                   filter_equal_dict={"status_axon":"non"})[["pt_root_id"]]
    client.version = curr_ver
prev_unpf_df = client.materialize.tables.nucleus_detection_v0() \
    .query(timestamp=client.materialize.get_version_metadata(prev_ver)["time_stamp"])[["pt_root_id"]]
prev_unpf_df = prev_unpf_df[~prev_unpf_df["pt_root_id"].isin(prev_proof_df["pt_root_id"])]

# update cell ids of previous version"s dataframe to latest ids
start = time.process_time()

try:
    with open(f"./results/id_update_dict_unpf_prev_ver={prev_ver}_num_updt={num_updt}.pkl", "rb") as f:
        update_dict = pickle.load(f)
except:
    prev_root_ids = prev_unpf_df["pt_root_id"].values
    expired_ids = prev_root_ids[~client.chunkedgraph.is_latest_roots(prev_root_ids)]
    
    rng = np.random.default_rng(0)
    expired_ids = rng.choice(expired_ids, size=min(num_updt, len(expired_ids)), replace=False)
    
    updated_ids = np.zeros_like(expired_ids)
    print(f"Updating {len(expired_ids)} cell ids to latest roots...")
    for i,root_id in tqdm(enumerate(expired_ids)):
        try:
            updated_ids[i] = client.chunkedgraph.suggest_latest_roots(root_id)
        except:
            updated_ids[i] = 0  # If the root ID is not found, set to 0
    update_dict = dict(zip(expired_ids, updated_ids))
    with open(f"./results/id_update_dict_unpf_prev_ver={prev_ver}_num_updt={num_updt}.pkl", "wb") as f:
        pickle.dump(update_dict, f)
prev_unpf_df["pt_root_id"] = prev_unpf_df["pt_root_id"].replace(update_dict)

print(f"Time taken to update previous unproofread dataframe: {time.process_time() - start} seconds")

# find shared cell ids between current and previous unproofread dataframes
shared_unpf_df = curr_unpf_df.merge(prev_unpf_df, on="pt_root_id", suffixes=("_curr","_prev"))

# add cell type and column information to shared dataframe
shared_unpf_df = shared_unpf_df.merge(type_df, on="pt_root_id", how="left")
shared_unpf_df["in_column"] = shared_unpf_df["pt_root_id"].isin(column_df["pt_root_id"])

# sample presynaptic cells, and compare synapses identities between current and previous versions
def get_synapses(pre_id):
    curr_out_syn_df = client.materialize.synapse_query(pre_ids=pre_id,remove_autapses=True)
    if client.chunkedgraph.is_latest_roots(pre_id,timestamp=client.materialize.get_timestamp(prev_ver)):
        prev_out_syn_df = client.materialize.synapse_query(pre_ids=pre_id,remove_autapses=True,
                                                           materialization_version=prev_ver)
    else:
        prev_out_syn_df = client.materialize.synapse_query(pre_ids=client.chunkedgraph.suggest_latest_roots(pre_id,
                                                           timestamp=client.materialize.get_timestamp(prev_ver)),
                                                           remove_autapses=True, materialization_version=prev_ver)
    expired_ids = np.unique(prev_out_syn_df["post_pt_root_id"]\
        .values[~client.chunkedgraph.is_latest_roots(prev_out_syn_df["post_pt_root_id"])])
    updated_ids = np.zeros_like(expired_ids)
    print(f"Updating {len(expired_ids)} post synaptic ids to latest roots...")
    for i,id in tqdm(enumerate(expired_ids)):
        try:
            updated_ids[i] = client.chunkedgraph.suggest_latest_roots(id)
        except:
            updated_ids[i] = 0  # If the root ID is not found, set to 0
    update_dict = dict(zip(expired_ids, updated_ids))
    prev_out_syn_df["post_pt_root_id"] = prev_out_syn_df["post_pt_root_id"].replace(update_dict)
    
    return curr_out_syn_df["post_pt_root_id"], prev_out_syn_df["post_pt_root_id"]

curr_num_syn = np.zeros(num_samp)
prev_num_syn = np.zeros(num_samp)
shared_syns = np.zeros(num_samp)

rng = np.random.default_rng(0)

start = time.process_time()

idxs = rng.choice(len(shared_unpf_df), size=num_samp, replace=False)
for i,idx in enumerate(idxs):
    print(f"Processing index {i+1}/{num_samp}...")
    curr_post_ids, prev_post_ids = get_synapses(shared_unpf_df["pt_root_id"].values[idx])
    curr_num_syn[i] = len(curr_post_ids)
    prev_num_syn[i] = len(prev_post_ids)
    shared_syns[i] = curr_post_ids.isin(prev_post_ids).sum()
    
print(f"Time taken to process synapses: {time.process_time() - start} seconds")
    
samp_unpf_df = shared_unpf_df.iloc[idxs]
samp_unpf_df["curr_num_syn"] = curr_num_syn
samp_unpf_df["prev_num_syn"] = prev_num_syn
samp_unpf_df["shared_syns"] = shared_syns

samp_unpf_df.to_csv(f"./results/unpf_syn_rel_prev_ver={prev_ver}_num_samp={num_samp}_seed={seed}.csv", index=False)

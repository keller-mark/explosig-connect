from explosig_connect import connect
import pandas as pd


conn = connect(session_id="4392af8f", empty=True, how=None, server_hostname='http://localhost:8100', client_hostname='http://localhost:8000')


# Load data
data_url = "https://raw.githubusercontent.com/keller-mark/explosig-connect/master/examples/data/{filename}"
# - counts data
sbs_counts_nosplit_df = pd.read_csv(data_url.format(filename="counts.DFCI-30-Kasar2015.nosplit.WGS.SBS-96.tsv"), sep='\t', index_col=0)
sbs_counts_nmd1000_df = pd.read_csv(data_url.format(filename="counts.DFCI-30-Kasar2015.nmd1000.WGS.SBS-96.tsv"), sep='\t', index_col=0)
# - signatures and exposures data,
#   preprocessed using the code here: https://github.com/keller-mark/Reproducing-Kasar2015
sbs_sigs_nosplit_df = pd.read_csv(data_url.format(filename="nosplit_run_5_W.txt"), sep='\t', index_col=0).transpose()
sbs_sigs_nmd1000_df = pd.read_csv(data_url.format(filename="nmd1000_run_30_W.txt"), sep='\t', index_col=0).transpose()
sbs_exps_nosplit_df = pd.read_csv(data_url.format(filename="nosplit_run_5_H.txt"), sep='\t', index_col=0).transpose()
sbs_exps_nmd1000_df = pd.read_csv(data_url.format(filename="nmd1000_run_30_H.txt"), sep='\t', index_col=0).transpose()


# Generate mutation type counts df
counts_nosplit_df = pd.DataFrame(index=sbs_counts_nosplit_df.index.values.tolist(), columns=["SBS", "DBS", "INDEL"], data=[{"SBS": sbs_count, "DBS": 0, "INDEL": 0} for sbs_count in sbs_counts_nosplit_df.sum(axis=1).values.tolist()])
counts_nmd1000_df = pd.DataFrame(index=sbs_counts_nmd1000_df.index.values.tolist(), columns=["SBS", "DBS", "INDEL"], data=[{"SBS": sbs_count, "DBS": 0, "INDEL": 0} for sbs_count in sbs_counts_nmd1000_df.sum(axis=1).values.tolist()])


# Generate samples metadata df
sample_metadata_nosplit_df = pd.DataFrame(
    index=sbs_counts_nosplit_df.index.values.tolist(),
    columns=["Study"],
    data=[ {"Study": "DFCI-30-Kasar2015"} for sample_id in sbs_counts_nosplit_df.index.values.tolist() ]
)
sample_metadata_nmd1000_df = pd.DataFrame(
    index=sbs_counts_nmd1000_df.index.values.tolist(),
    columns=["Study"],
    data=[ {"Study": "DFCI-30-Kasar2015"} for sample_id in sbs_counts_nmd1000_df.index.values.tolist() ]
)

def send_nosplit_data():
    conn.send_sample_metadata(sample_metadata_nosplit_df)
    conn.send_mutation_type_counts(counts_nosplit_df)
    conn.send_signatures('SBS', sbs_sigs_nosplit_df)
    conn.send_exposures('SBS', sbs_exps_nosplit_df)

def send_nmd1000_data():
    conn.send_sample_metadata(sample_metadata_nmd1000_df)
    conn.send_mutation_type_counts(counts_nmd1000_df)
    conn.send_signatures('SBS', sbs_sigs_nmd1000_df)
    conn.send_exposures('SBS', sbs_exps_nmd1000_df)


# Send non-split data
#send_nosplit_data()

# Send split data
send_nmd1000_data()


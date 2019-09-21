from explosig_connect import connect
import pandas as pd

conn = connect(session_id="4392af8f", empty=True, how=None, server_hostname='http://localhost:8100', client_hostname='http://localhost:8000')

df = pd.DataFrame(index=["SA 01", "SA 02", "SA 03"], columns=["Study", "Donor"], data=[
    {"Study": "Breast", "Donor": "DO 01"},
    {"Study": "Breast", "Donor": "DO 02"},
    {"Study": "Pancreatic", "Donor": "DO 02"}
])
conn.send_sample_metadata(df)

df = pd.DataFrame(index=["SA 01", "SA 02", "SA 03"], columns=["SBS", "DBS", "INDEL"], data=[
    {"SBS": 10, "DBS": 10, "INDEL": 10},
    {"SBS": 10, "DBS": 10, "INDEL": 1},
    {"SBS": 10, "DBS": 10, "INDEL": 1}
])
conn.send_mutation_type_counts(df)
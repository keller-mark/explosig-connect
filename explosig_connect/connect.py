import requests
from .connection import EmptyConnection, ConfigConnection

HOSTNAME = 'https://explosig.lrgr.io'

def login(password, hostname):
    r = requests.post(hostname + '/login')
    r.raise_for_status()
    return r.json()['token']

def connect(session_id=None, password=None, hostname=HOSTNAME, how='auto'):
    if password != None and hostname != HOSTNAME:
        token = login(password, hostname)
    else:
        token = None

    if session_id == None:
        conn = EmptyConnection(token, hostname)
        conn.open(how=how)
        return conn
    else:
        return ConfigConnection(session_id, token, hostname)

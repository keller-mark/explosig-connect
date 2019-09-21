import requests
from .connection import EmptyConnection, ConfigConnection

CLIENT_HOSTNAME = 'https://explosig.lrgr.io'
SERVER_HOSTNAME = 'https://explosig-server.lrgr.io'

def login(password, server_hostname):
    r = requests.post(server_hostname + '/login')
    r.raise_for_status()
    return r.json()['token']

def connect(session_id=None, empty=False, password=None, server_hostname=SERVER_HOSTNAME, client_hostname=CLIENT_HOSTNAME, how='auto'):
    if password != None and server_hostname != SERVER_HOSTNAME:
        token = login(password, server_hostname)
    else:
        token = None

    if session_id == None or empty == True:
        conn = EmptyConnection(session_id, token, server_hostname, client_hostname)
        conn.open(how=how)
        return conn
    else:
        return ConfigConnection(session_id, token, server_hostname, client_hostname)

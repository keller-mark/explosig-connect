import requests
from .connection import EmptyConnection, ConfigConnection

CLIENT_HOSTNAME = 'https://explosig.lrgr.io'
SERVER_HOSTNAME = 'https://explosig-server.lrgr.io'

def _login(password, server_hostname):
    r = requests.post(server_hostname + '/login')
    r.raise_for_status()
    return r.json()['token']

def connect(session_id=None, empty=False, password=None, server_hostname=SERVER_HOSTNAME, client_hostname=CLIENT_HOSTNAME, how='auto'):
    """Connect to an ExploSig session.
    
    Parameters
    ----------
    session_id : `str`, optional
        An ExploSig session ID.
        If not provided, a new "empty" session will be started. 
        by default `None`
    empty : `bool`, optional
        If `True`, will open an "empty" session regardless of whether 
        `session_id` is provided.
        by default `False`
    password : `str`, optional
        Required if not using the public instance of ExploSig and instead using
        a password-protected instance.
        by default `None`
    server_hostname : `str`, optional
        Use to specify an alternate ExploSig server instance.
        by default `'https://explosig-server.lrgr.io'`
    client_hostname : `str`, optional
        Use to specify an alternate ExploSig client instance.
        by default `'https://explosig.lrgr.io'`
    how : `str`, optional
        If starting a new empty session, 
        the method for opening in the browser,
        passed to `explosig_connect.connection.EmptyConnection.open`.
        by default `'auto'`
    
    Returns
    -------
    `explosig_connect.connection.Connection`
        Returns an object of a Connection subclass: `explosig_connect.connection.EmptyConnection` 
        if starting a new "empty" session, or `explosig_connect.connection.ConfigConnection` 
        if connecting to an existing session that has been configured and started 
        from within ExploSig.
    """
    if password != None and server_hostname != SERVER_HOSTNAME:
        token = _login(password, server_hostname)
    else:
        token = None

    if session_id == None or empty == True:
        conn = EmptyConnection(session_id, token, server_hostname, client_hostname)
        conn.open(how=how)
        return conn
    else:
        return ConfigConnection(session_id, token, server_hostname, client_hostname)

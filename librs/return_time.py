from datetime import datetime as dt


def return_time(message: str = None):
    t = dt.now()
    s = str(t).split('.')[0]
    if message is None:
        print("[%s] Return time: no message to print" % s)
    else:
        print("[%s] %s" % (s, message))

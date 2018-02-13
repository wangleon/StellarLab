
def _str_to_float(string, exception_value=None):
    '''
    Convert string to float. Return `exception_value` if failed.
    '''
    try:
        return float(string)
    except:
        return exception_value

def _str_to_int(string, exception_value=None):
    '''
    Convert string to integer. Return `exception_value` if failed.
    '''
    try:
        return int(string)
    except:
        return exception_value

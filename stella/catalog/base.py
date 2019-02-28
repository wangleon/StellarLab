
def _str_to_float(string, exception=None):
    """Convert string to float. Return `exception_value` if failed.

    Args:
        string (str): Input string.
        exception (): Exceptional value returned if failed.

    Returns:
        int:
    """
    try:
        return float(string)
    except:
        return exception

def _str_to_int(string, exception=None):
    """Convert string to integer. Return `exception_value` if failed.

    Args:
        string (str): Input string.
        exception (): Exceptional value returned if failed.

    Returns:
        int:
    """
    try:
        return int(string)
    except:
        return exception

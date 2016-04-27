class FileNotFound(Exception):
    def __init__(self, filename):
        self.filename = filename
    def __str__(self):
        return repr(self.filename)

class ItemNotFound(Exception):
    def __init__(self, itemid):
        self.itemid = itemid
    def __str__(self):
        return repr(self.itemid)

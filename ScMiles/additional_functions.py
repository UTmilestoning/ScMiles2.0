def create_folder(path):
    import os
    if not os.path.exists(path):
        os.makedirs(path)

def ms_to_path(ms):
    import re
    lst = list(map(int, (re.findall('\d+', ms))))
    return lst

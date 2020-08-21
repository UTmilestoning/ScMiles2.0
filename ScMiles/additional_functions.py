def create_folder(path):
    import os
    if not os.path.exists(path):
        os.makedirs(path)


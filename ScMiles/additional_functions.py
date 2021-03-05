def create_folder(path):
    '''
    If a folder does not exist, the folder is created. If it already exists,
    it stays as is
    '''
    import os
    if not os.path.exists(path):
        os.makedirs(path)

def get_anchors(anchors):
    '''
    This changes a milestone in the form of "MS1_2" into anchor form,
    so it returns a list of the anchors (with each anchor as an integer), such as [1,2]
    '''
    import re
    if type(anchors) == int:
        return [anchors, None]
    else:
        return list(map(int, (re.findall('\d+', anchors))))

def get_initial_ms(path):
    '''get initial milestone for each free trajectory based on the folder name'''
    import re
    path_split = path.split("/")
    initial_ms = list(map(int,(re.findall('\d+', path_split[-3]))))
    with open(path + '/start.txt', 'w+') as f1:
        print(initial_ms[0], initial_ms[1], file=f1)    
    return initial_ms

def get_next_frame_num(struPath):
    import os
    '''
    Description: This function is used to find the next non-existent folder
        So if we have the folders 
        crd/1_2/1, crd/1_2/2, crd/1_2/3
        This would return 4 since it is the first one to not exist.
    Arguments: struPath, the path that we are looking at. So for the example in the 
        description, this would be 'crd/1_2'
    Returns: next_frame, which is the next non-existent folder (4 in the example)
    '''
    next_frame = 1
    while True:
        pdbPath = struPath + '/' + str(next_frame) 
        if os.path.exists(pdbPath):
            next_frame += 1
        else:
            return next_frame

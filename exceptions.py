class ImproperStructureFiles(Exception):
    '''The given structure files are not useable.'''
    pass

class ImproperInputType(Exception):
    '''The given input type is invalid.'''
    pass

class ImproperInputRange(Exception):
    '''The range of inputs is invalid.'''
    pass

class ImproperMutations(Exception):
    '''The mutations file is invalid.'''
    pass

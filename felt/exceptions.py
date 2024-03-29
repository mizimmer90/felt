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

class InvalidData(Exception):
    '''The provided data is not valid.'''
    pass

class PathExists(Exception):
    '''The specified path already exists.'''
    pass

class ImproperDihedralRotation(Exception):
    '''The specified dihedral cannot be rotated.'''
    pass

class UnexpectedError(Exception):
    '''There was an unexpected error'''
    pass

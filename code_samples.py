from realize import Realizable_Set

# Here's how can check if a set is realizable:
K = Realizable_Set.get_instance((2,3))
x=K.is_realizable() # should return `True`

K = Realizable_Set.get_instance((1,2))
x=K.is_realizable() # should return `False`

# You can check if a set is blocking:
K = Realizable_Set.get_instance((1,))
x=K.is_blocking() # should return `True`

K = Realizable_Set.get_instance((2,3))
x=K.is_blocking() # should return `False`

K = Realizable_Set.get_instance((2,3,4))
x=K.is_blocking() # should return `True`

# You can calculate every single permutation matrix that realizes a given set.
# There are two methods that do that, `get_realizations` and `get_blocks`. The first returns p x p matrices,
# where p is the maximal element of K, and the second returns n x n matrices, where n is the number of
# elements of K.

K = Realizable_Set.get_instance((2,3,4))
x=list(K.get_realizations()) # returns a `list` of 5 x 5 matrices. There should be exactly 4 matrices.
x=K.get_blocks() # returns a `list` of 3 x 3 matrices. There should be exactly 2 matrices.

K = Realizable_Set.get_instance((3,4,5,6,7))
list(K.get_realizations()) # returns a `list` of 8 x 8 matrices. There should be exactly 36 matrices.
K.get_blocks() # returns a `list` of 5 x 5 matrices. There should be exactly 6 matrices.

K = Realizable_Set.get_instance((1,4))
list(K.get_realizations()) # returns a `list` of 5 x 5 matrices. There should be exactly 8 matrices.

K = Realizable_Set.get_instance((3,4,5,6))
list(K.get_realizations()) # returns a `list` of 7 x 7 matrices. There should be exactly 60 matrices.

K = Realizable_Set.get_instance((3,))
list(K.get_realizations()) # returns a `list` of 4 x 4 matrices. There should be exactly 3 matrices.

K = Realizable_Set.get_instance((4,))
list(K.get_realizations()) # returns a `list` of 5 x 5 matrices. There should be exactly 4 matrices.

# You can calculate all the *symmetric* matrices realizing a given set.
# There is only one method that does that, `get_symmetric_realizations`. This method takes an additional argument
# specifying the size of symmetric matrix.
K = Realizable_Set.get_instance((4,5,6))
K.get_symmetric_realizations(7) # returns a `list` of 7 x 7 matrices. There should be exactly 6 matrices.





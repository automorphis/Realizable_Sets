from main import Realizable_Set

# Here's how can check if a set is realizable:
K = Realizable_Set.get_instance((2,3))
K.is_realizable() # should return `True`

K = Realizable_Set.get_instance((1,2))
K.is_realizable() # should return `False`

# You can check if a set is blocking:
K = Realizable_Set.get_instance((1,))
K.is_blocking() # should return `True`

K = Realizable_Set.get_instance((2,3))
K.is_blocking() # should return `False`

K = Realizable_Set.get_instance((2,3,4))
K.is_blocking() # should return `True`

# You can calculate every single permutation matrix that realizes a given set.
# There are two methods that do that, `get_realizations` and `get_blocks`. The first returns p x p matrices, where p
# is the maximal element of K, and the second returns n x n matrices, where n is the number of elements of K.
#
# The method `get_realizations` is slightly broken at this moment. It does correctly calculate realizations, but I am
# not happy with the form of the output. I'll try to fix it soon.
K = Realizable_Set.get_instance((2,3,4))
K.get_realizations() # returns a `list` of 4 x 4 matrices. There should be exactly 2 matrices.
K.get_blocks() # returns a `list` of 3 x 3 matrices. There should be exactly 2 matrices.

K = Realizable_Set.get_instance((3,4,5,6,7))
K.get_realizations() # returns a `list` of 7 x 7 matrices. There should be exactly 6 matrices.
K.get_blocks() # returns a `list` of 5 x 5 matrices. There should be exactly 6 matrices.

# You can calculate all the *symmetric* matrices realizing a given set.
# There is only one method that does that, `get_symmetric_realizations`. This method takes an additional argument
# specifying the size of symmetric matrix.
K = Realizable_Set.get_instance((4,5,6))
K.get_symmetric_realizations(7) # returns a `list` of 7 x 7 matrices. There should be exactly 6 matrices.





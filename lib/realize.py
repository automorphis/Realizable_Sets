import copy

import numpy as np
from itertools import chain, combinations, product, accumulate


class Not_Realizable (RuntimeError):
    """Raised if a given set is not realizable."""
    pass

class Not_Blocking (RuntimeError):
    """Raised if a given set is not blocking."""
    pass

class Not_Realizing_Matrix (RuntimeError):
    """Raised if a given matrix does not realize any set."""
    pass

class Not_Realizing_Matrix_Template (RuntimeError):
    """Raised if a given matrix is not a realizing template for any set."""
    pass

class Not_Instantiated_Error (RuntimeError):
    """Raised if the methods `set_array` or `set_cols` are not called before calling any other methods of
    `Realizing_Matrix`."""
    def __init__(self, method_name):
        super().__init__(
            f"You must call either the method `set_cols` or `set_array` before calling the method `{method_name}`."
        )

class Already_Instantiated_Error (RuntimeError):
    """Raised if either of the methods `set_array` or `set_cols` are called twice on a single `Realizing_Matrix`
    instance."""
    def __init__(self):
        super().__init__(
            "You cannot call either of the methods `set_array` or `set_cols` twice on a single `Realizing_Matrix` " +
            "instance."
        )

class Realizable_Set:
    """This class represents a realizable set.

    This class memoizes its instances, as long as `Realizable_Set.remember_instances is True` (the default value is
    `True`). Memoized instances can be accessed via the static method `get_instance`.
    """

    instances = {}

    remember_instances = True

    @staticmethod
    def get_instance(K):
        """Return a memoized `Realizable_Set` instance.

        Successfully calling this method (or the `Realizable_Set` constructor method) is not a guarantee that `K`
        represents a `Realizable_Set`. Call the method `is_realizable` to check if your instance is indeed realizable.

        :param K: An iterable of positive integers.
        :return: type `Realizable_Set`.
        """

        if not Realizable_Set.remember_instances:
            return Realizable_Set(K)

        else:
            K = tuple(sorted(list(K)))
            if K in Realizable_Set.instances.keys():
                return Realizable_Set.instances[K]
            else:
                ret = Realizable_Set(K)
                Realizable_Set.instances[K] = ret
                return ret

    def __init__(self, K):
        """Using this constructor method bypasses memoization. It is preferred to use the `Realizable_Set.get_instance`
        static method, which will call this constructor method if necessary.

        Successfully calling this method (or the `Realizable_Set.get_instance` static method) is not a guarantee that
        `K` represents a `Realizable_Set`. Call the method `is_realizable` to check if your instance is indeed
        realizable.

        :param K: An iterable of positive `int`s.
        """

        K = list(K)
        for k in K:
            if not isinstance(k, (int, np.int32, np.int64)):
                raise TypeError("elements of `K` must be of type `int`")
            elif k < 0:
                raise ValueError("elements of `K` must be positive")

        self.K = tuple(sorted(K))

        if len(set(self.K)) != len(self.K):
            raise ValueError("elements of `K` must be unique.")

        self._is_realizable_nece_cond = None
        self._blocks = None
        self._num_blocks = None
        self._realization_templates = None
        self._symm_realizing_mats = None
        self._realizable_suff_cond = None
        self._is_realizable = None
        self._is_blocking = None
        self._skinny_dim = None
        self._fat_dim = None
        self._large_blking_pfx = None
        self._norm_compl_large_blking_pfx = None

    def __len__(self):
        return len(self.K)

    def __iter__(self):
        return iter(self.K)

    def __getitem__(self, index):
        """Enumerate K as k_0 < k_1 < ... k_{n-1}. This method returns k_{index}.

        :param index: Non-negative `int`.
        :return: Positive `int`.
        """
        return self.K[index]

    def __contains__(self, item):
        return item in self.K

    def __hash__(self):
        return hash(self.K)

    def __eq__(self, other):
        return self.K == other.K

    def __repr__(self):
        return str(self.K)

    def __str__(self):
        return repr(self)

    def __format__(self, format_spec):
        """I don't remember if this method works properly and I don't care to check in this moment."""
        return self.K.__format__(format_spec)

    def get_max(self):
        return self[-1]

    def get_min(self):
        return self[0]

    def index(self, item):
        return self.K.index(item)

    def get_initial_segment(self, i):
        """Return {k_1, ..., k_i}, where `n = len(self` and k_1 < ... < k_i < ... < k_n are the ordered elements
        of this instance.

        :param i: (type `int`) positive integer.
        :raise TypeError: If `i` is not type `int`.
        :raise ValueError: If `i` is too large or too small.
        :return: (type `Realizable_Set`)
        """

        if not isinstance(i, int):
            raise TypeError("`i` is not an `int`")

        if not(1 <= i <= len(self)):
            raise ValueError("`i` must be between 1 and `len(self)`")

        return Realizable_Set.get_instance(self.K[:i])

    def is_realizable_necessary_condition(self):
        """Check if this instance satisfies the realizability necessary condition, as follows: Let `n = len(K)` and
        write K = {k_1, ..., k_n}, where k_1 < ... < k_n. For all i = 1, ..., n-1, if {k_1, ..., k_i} is blocking, then
        K \cap [k_i + 1, 2n] is empty.

        The necessary condition above is not sufficient, so non-realizable sets may return `True`.

        It is best to call this function when `Realizable_Set.remember_instances is True`, otherwise it can be quite
        slow.

        :return: (type `bool`)

        """

        if self._is_realizable:
            self._set_is_realizable_nece_cond(True)
            return True

        elif self._is_realizable_nece_cond is not None:
            return self._is_realizable_nece_cond

        else:

            for i, k in zip(range(1, len(self)), self.K[:-1]):
                Kp = self.get_initial_segment(i)
                if Kp.is_blocking(False) and any(l in range(k+1,2*i+1) for l in self.K[i:]):
                    self._set_is_realizable_nece_cond(False)
                    break

            else:
                self._set_is_realizable_nece_cond(True)

            return self._is_realizable_nece_cond

    def is_realizable(self, use_necessary_condition = True):
        """Check if this instance is indeed realizable.

        :param use_necessary_condition: (type `bool`, default `True`). Whether to call the method
        `is_realizable_necessary_condition` to check for non-realizability. Warning: That method can be quite slow if
        `Realizable_Set.remember_instances is False`, so it is best to set `use_necessary_condition = False` in that
        case.
        :return: (type `bool)

        """

        if self._is_realizable is not None:
            return self._is_realizable

        else:

            if len(self) == 1:
                self._set_is_realizable(True)
                return self._is_realizable

            if Realizable_Set.remember_instances:
                realizable = True
                for i in range(len(self)):
                    Kp = self.get_initial_segment(i+1)

                    if not realizable:
                        Kp._set_is_realizable(False)

                    elif Kp._is_realizable is None:

                            nece_cond = True
                            if use_necessary_condition:
                                nece_cond = Kp.is_realizable_necessary_condition()

                            if nece_cond:

                                Ki = self[i]
                                Kp._set_is_realizable(Kp._is_realizable_dfs_helper(
                                    0,
                                    [None] * (Ki + 1),
                                    list(range(Ki + 1 - self[0], Ki + 1))
                                ))
                                if not Kp._is_realizable:
                                    realizable = False

                            else:
                                Kp._set_is_realizable(False)
                                realizable = False

                    elif not Kp._is_realizable:
                        realizable = False

                return self._is_realizable

            else:

                nece_cond = True
                if use_necessary_condition:
                    nece_cond = self.is_realizable_necessary_condition()

                if nece_cond:

                    kp = self.get_max()
                    self._set_is_realizable(self._is_realizable_dfs_helper(
                        0,
                        [None] * (kp+1),
                        list(range(kp+1 - self[0], kp+1))
                    ))

                else:
                    self._set_is_realizable(False)

                return self._is_realizable

    def _set_is_realizable(self, value):
        self._is_realizable = value
        if value:
            self._is_realizable_nece_cond = True
        else:
            self._set_num_blocks(0)

    def _set_is_realizable_nece_cond(self, value):
        self._is_realizable_nece_cond = value
        if not value:
            self._is_realizable = False
            self._set_num_blocks(0)

    def get_conjectured_realizable_necessary_condition_vector(self):
        """If the realizable set K is enumerated as k_0 < k_1 < ... < k_{n-1}, then this method returns the list of
        differences k_i - (3/2)i + 1 for all i = 0, 1, ..., n-1. It is conjectured (Conjecture 0.6) that this list is
        non-negative.

        :return: A `list` of `int`s.
        """
        return [k - ((3 / 2) * i + 1) for i,k in enumerate(self)]

    def get_realization_necessary_condition_direction_vector(self):
        vect = self.get_conjectured_realizable_necessary_condition_vector()
        return ["u" if v0 <= v1 else "d" for v0,v1 in zip(vect[:-1],vect[1:])]

    def get_realization_templates(self):
        """Calculate all (k+1) x (k+1) strictly upper triangular matrices that realize this `Realizable_Set`, where k
        is the largest element of this instance.

        If `Realizing_Matrix_Template.remember_instances is True`, then this instance will remember the output of this
        method and return the same `list` on future calls. Otherwise, this method will recalculate its output.

        :return: A `list` of `Realizing_Matrix_Template`. If no matrices realize this set, an empty list is returned.
        """

        kp = self.get_max()
        size = kp + 1

        if self._realization_templates is not None:
            return self._realization_templates

        else:

            try:
                ret = self._get_realization_templates_dfs_helper(
                    0,
                    [None]*size,
                    list(range(size - self.get_min(), size))
                )

            except Not_Realizable:
                ret = []

            if Realizing_Matrix_Template.remember_instances:
                self._realization_templates = ret

            return ret

    def get_realizations(self):
        """Yield all (k+1) x (k+1) permutation matrices that realize this instance, assuming they exist, where k is the
        maximum element of this `Realizable_Set`.

        This method does NOT remember its output, even if `Realizing_Matrix.remember_instances is True`, due to the
        huge amount of memory that would require.

        Warning: If k is large relative to n, where n is `len(self)`, then this generator may return a huge number of
        elements.
        """
        for template in self.get_realization_templates():
            for real in template.iterate_realizations():
                yield real

    def is_blocking(self, use_necessary_condition = True, use_conjectured_shortcircuit = False):
        """Check if this instance represents a blocking set.

        :param use_necessary_condition: (type `bool`, default `True`). Whether to call the method
        `is_realizable_necessary_condition` to check for non-realizability. Warning: That method can be quite slow if
        `Realizable_Set.remember_instances is False`, so it is best to set `use_necessary_condition = False` in that
        case.
        :param use_conjectured_shortcircuit: (type `bool`, default `False`). Return `True` if there at least one
        blocking realization, rather than check that all realizations are blocking.
        :return: `True` if this instance represents a blocking set, and `False` otherwise.
        """

        if self._is_blocking is not None:
            return self._is_blocking

        elif use_conjectured_shortcircuit:
            kp = self.get_max()
            size = kp + 1
            try:
                self._set_is_blocking(self._is_blocking_shortcircuit_dfs_helper(
                    0,
                    [None] * size,
                    list(range(size - self.get_min(), size))
                ))
            except Not_Realizable:
                self._set_is_realizable(False)

            return self._is_blocking

        elif use_necessary_condition:

            if not self.is_realizable_necessary_condition():
                self._set_is_blocking(False)
                if Realizing_Matrix.remember_instances:
                    self._blocks = []
                return False

            _, K = self.get_largest_blocking_prefix()

            if K is not None:
                self._set_is_blocking(K.is_blocking(False))
            else:
                self._set_is_blocking(len(self.get_blocks(False)) > 0)

            return self._is_blocking

        else:
            self._set_is_blocking(len(self.get_blocks(False)) > 0)
            return self._is_blocking

    def get_blocks(self, use_necessary_condition = True):
        """If this instance represents a blocking set, return all the matrices that realize it. If not, then
        return an empty `list`.

        If `Realizing_Matrix.remember_instances is True`, then this instance will remember the output of this
        method and return the same `list` on future calls. Otherwise, this method will recalculate its output.

        :param use_necessary_condition: (type `bool`, default `True`). Whether to call the method
        `is_realizable_necessary_condition` to check for non-realizability. Warning: That method can be quite slow if
        `Realizable_Set.remember_instances is False`, so it is best to set `use_necessary_condition = False` in that
        case.
        :return: A `list` of `Realizing_Matrix`.
        """

        if self._blocks is not None:
            return self._blocks

        elif self._is_realizable is not None and not self._is_realizable:
            ret = []
            if Realizing_Matrix.remember_instances:
                self._blocks = ret
            self._set_num_blocks(0)
            return ret

        else:

            n = len(self)
            kp = self.get_max()
            size = kp + 1

            if 2 * n <= kp:
                ret = []

            else:

                if use_necessary_condition and not self.is_realizable_necessary_condition():
                    ret = []
                    is_blocking = False

                else:

                    try:
                        ret = self._get_blocks_dfs_helper(
                            0,
                            [None] * size,
                            list(range(size - self.get_min(), size))
                        )
                        is_blocking = True

                    except (Not_Realizable, Not_Blocking) as e:
                        if isinstance(e, Not_Realizable):
                            self._set_is_realizable(False)
                        ret = []
                        is_blocking = False

                if is_blocking:
                    blocks = []
                    for cols in ret:
                        mat = Realizing_Matrix()
                        mat.set_cols([c + n - size for c in cols[:n]], kp, True)
                        mat = Realizing_Matrix.get_instance(mat)
                        blocks.append(mat)
                    self._set_is_realizable(True)
                    ret = blocks

            if Realizing_Matrix.remember_instances:
                self._blocks = ret

            self._set_num_blocks(len(ret))

            return ret

    def get_num_blocks(self, use_necessary_condition = True):

        if self._num_blocks is not None:
            return self._num_blocks

        elif self._blocks is not None:
            self._num_blocks = len(self._blocks)
            return self._num_blocks

        elif self._is_realizable is not None and not self._is_realizable:
            self._num_blocks = 0
            return 0

        else:

            if use_necessary_condition:

                if not self.is_realizable_necessary_condition():
                    self._set_num_blocks(0)
                    if Realizing_Matrix.remember_instances:
                        self._blocks = []
                    return 0

                K, Kc = self.get_largest_blocking_prefix()

                if K is not None:
                    self._set_num_blocks(K.get_num_blocks(False) * Kc.get_num_blocks(True))

                else:
                    self._set_num_blocks(len(self.get_blocks()))

            else:

                self._set_num_blocks(len(self.get_blocks()))

            return self._num_blocks

    def get_largest_blocking_prefix(self):
        """This method returns the largest contiguous subset {k_1, ..., k_i} that is blocking, as well as its
        normalized complement.

        Even if this instance is not realizable, this method will still return.

        :return: Tuple of length 2, first element is the largest blocking prefix, second is the "normalized complement"
        of the largest blocking prefix, defined as {k_{i+1} - 2i, ..., k_n - 2i}. If no prefix is blocking, returns
        (None, None).
        """

        if self._large_blking_pfx is not None:
            if self._large_blking_pfx != 0:
                return self._large_blking_pfx, self._norm_compl_large_blking_pfx
            else:
                return None, None

        else:

            n = len(self)
            for i in range(n - 1, 0, -1):
                Kp = self.get_initial_segment(i)
                if Kp.is_blocking(False):
                    K = Realizable_Set.get_instance([k - 2 * i for k in self.K[i:]])
                    break
            else:
                Kp = K = None

            if Kp is not None:
                self._large_blking_pfx = Kp
                self._norm_compl_large_blking_pfx = K
            else:
                self._large_blking_pfx = self._norm_compl_large_blking_pfx = 0

            return Kp, K

    def _set_num_blocks(self, num):
        self._num_blocks = num
        self._is_blocking = num > 0

    def _set_is_blocking(self, value):
        self._is_blocking = value
        if not value:
            self._num_blocks = 0

    def get_symmetric_realizations(self, m):
        """Calculate all m x m symmetric permutation matrices that realize this instance.

        :param m: Positive `int`.
        :return: A `list` of `Realizing_Matrix`.
        """

        if self._symm_realizing_mats is not None:
            return self._symm_realizing_mats

        else:

            try:
                realizing_cols = self._get_symmetric_realizations_dfs_helper(
                    m,
                    0,
                    [None]*m,
                    list(range(m-self[0], m))
                )
            except Not_Realizable:
                pass

            self._symm_realizing_mats = []
            for cols in realizing_cols:
                for i,j in enumerate(cols):
                    if j is not None:
                        cols[j] = i
                cols = [j if j is not None else i for i,j in enumerate(cols)]
                A = Realizing_Matrix()
                A.set_cols(cols,m-1,True)
                A = Realizing_Matrix.get_instance(A)
                self._symm_realizing_mats.append(A)

            return self._symm_realizing_mats

    def get_skinny_fat_rectangle_dimensions(self):

        if self._skinny_dim is not None:
            return self._skinny_dim, self._fat_dim

        else:

            real_temps = self.get_realization_templates()

            if not self.is_realizable():
                raise Not_Realizable

            kp = self.get_max()

            skinny_diff = 0
            fat_diff = kp - 1

            skinny_dim = (None, None)
            fat_dim = (None, None)

            for real_temp in real_temps:

                dim = real_temp.get_dims()
                diff = abs(dim[0] - dim[1])

                if (diff > skinny_diff or (
                    diff == skinny_diff and (
                            skinny_dim[0] is None or
                            min(*dim) > skinny_dim[0]
                        )
                    )
                ):
                    skinny_dim = (min(dim), max(dim))
                    skinny_diff = diff

                if (diff < fat_diff or (
                    diff == fat_diff and (
                            fat_dim[0] is None or
                            min(*dim) < fat_dim[0]
                        )
                    )
                ):
                    fat_dim = (min(dim), max(dim))
                    fat_diff = diff

            self._skinny_dim = skinny_dim
            self._fat_dim = fat_dim

            return self._skinny_dim, self._fat_dim

    def _get_pre_next_frontier_cols_helper(self, curr_k_index, cols):

        n = len(self)
        size = self.get_max() + 1

        if curr_k_index + 1 < n:
            next_k = self[curr_k_index + 1]
            return [
                j
                for j in range(size - next_k, size)
                if j not in cols and cols[j + next_k - size] is None
            ]

        else:
            return []

    def _get_next_cols_helper(self, curr_k_index, j, cols):

        size = self.get_max() + 1
        i = j - size + self[curr_k_index]
        next_cols = copy.copy(cols)
        next_cols[i] = j
        return next_cols

    def _get_next_frontier_cols_helper(self, curr_k_index, j, pre_next_frontier_cols):

        n = len(self)
        i = j - n + self[curr_k_index]

        if curr_k_index + 1 < n:
            next_frontier_cols = copy.copy(pre_next_frontier_cols)

            for jp in [j, n - self[curr_k_index + 1] + i]:
                try:
                    next_frontier_cols.remove(jp)
                except ValueError:
                    pass

            return next_frontier_cols

        else:
            return []

    def _is_realizable_dfs_helper(self, curr_k_index, cols, frontier_cols):

        n = len(self)

        if len(frontier_cols) == 0:
            return sum(j is not None for j in cols) == n

        pre_next_frontier_cols = self._get_pre_next_frontier_cols_helper(curr_k_index, cols)

        for j in frontier_cols:

            next_cols = self._get_next_cols_helper(curr_k_index, j, cols)

            next_frontier_cols = self._get_next_frontier_cols_helper(curr_k_index, j, pre_next_frontier_cols)

            if self._is_realizable_dfs_helper(
                curr_k_index + 1,
                next_cols,
                next_frontier_cols
            ):
                return True

        return False

    def _is_blocking_shortcircuit_dfs_helper(self, curr_k_index, cols, frontier_cols):

        kp = self.get_max()
        size = kp + 1
        n = len(self)

        if len(frontier_cols) == 0:

            if sum(j is not None for j in cols) == n:
                return None in cols and cols.index(None) == n and min(cols[:n]) == size - n
            else:
                raise Not_Realizable

        pre_next_frontier_cols = self._get_pre_next_frontier_cols_helper(curr_k_index, cols)

        non_realizable_branch = True

        for j in frontier_cols:

            next_cols = self._get_next_cols_helper(curr_k_index, j, cols)

            next_frontier_cols = self._get_next_frontier_cols_helper(curr_k_index, j, pre_next_frontier_cols)

            try:
                return self._is_blocking_shortcircuit_dfs_helper(
                    curr_k_index + 1,
                    next_cols,
                    next_frontier_cols
                )

            except Not_Realizable:
                pass

        raise Not_Realizable

    def _get_blocks_dfs_helper(self, curr_k_index, cols, frontier_cols):

        kp = self.get_max()
        size = kp + 1
        n = len(self)

        if len(frontier_cols) == 0:

            if sum(j is not None for j in cols) == n:

                if None in cols and cols.index(None) == n and min(cols[:n]) == size - n:
                    return [cols]
                else:
                    raise Not_Blocking

            else:
                raise Not_Realizable

        pre_next_frontier_cols = self._get_pre_next_frontier_cols_helper(curr_k_index, cols)

        non_realizable_branch = True
        realizing_cols = []

        for j in frontier_cols:

            next_cols = self._get_next_cols_helper(curr_k_index, j, cols)

            next_frontier_cols = self._get_next_frontier_cols_helper(curr_k_index, j, pre_next_frontier_cols)

            try:
                _realizing_cols = self._get_blocks_dfs_helper(
                    curr_k_index + 1,
                    next_cols,
                    next_frontier_cols
                )

                if len(_realizing_cols) > 0:
                    non_realizable_branch = False
                    realizing_cols.extend(_realizing_cols)

                else:
                    return []

            except Not_Realizable:
                pass

        if non_realizable_branch:
            raise Not_Realizable
        else:
            return realizing_cols

    def _get_symmetric_realizations_dfs_helper(self, m, curr_k_index, cols, frontier_cols):

        if len(frontier_cols) == 0:
            if sum(cols[i] is not None for i in range(m)) == len(self):
                return [cols]
            else:
                raise Not_Realizable

        curr_k = self[curr_k_index]
        next_k_index = curr_k_index + 1

        if next_k_index < len(self):
            next_k = self[next_k_index]
            pre_next_frontier_cols = [
                j
                for j in range( max(0, m-next_k),  min(m, 2*m - next_k))
                if j not in cols and j+next_k-m not in cols and cols[j + next_k - m] is None and cols[j] is None
            ]

        non_realizable_branch = True
        realizing_cols = []

        for j in frontier_cols:

            i = j - m + curr_k
            next_cols = copy.copy(cols)
            next_cols[i] = j

            if next_k_index < len(self):
                next_frontier_cols = copy.copy(pre_next_frontier_cols)
                for jp in [j, i, m-next_k+j, m-next_k+i]:
                    try:
                        next_frontier_cols.remove(jp)
                    except ValueError:
                        pass
            else:
                next_frontier_cols = []

            try:
                _realizing_cols = self._get_symmetric_realizations_dfs_helper(
                    m,
                    next_k_index,
                    next_cols,
                    next_frontier_cols
                )

                if len(_realizing_cols) > 0:
                    non_realizable_branch = False
                    realizing_cols.extend(_realizing_cols)

                else:
                    return []

            except Not_Realizable:
                pass

        if non_realizable_branch:
            raise Not_Realizable
        else:
            return realizing_cols

    def _get_realization_templates_dfs_helper(self, curr_k_index, cols, frontier_cols):

        kp = self.get_max()
        size = kp + 1

        if len(frontier_cols) == 0:

            if sum(j is not None for j in cols) == len(self):
                template = Realizing_Matrix_Template()
                template.set_cols(cols, kp, True)
                template = Realizing_Matrix_Template.get_instance(template)
                return [template]

            else:
                raise Not_Realizable

        pre_next_frontier_cols = self._get_pre_next_frontier_cols_helper(curr_k_index, cols)

        non_realizable_branch = True
        realizing_cols = []

        for j in frontier_cols:

            next_cols = self._get_next_cols_helper(curr_k_index, j, cols)

            next_frontier_cols = self._get_next_frontier_cols_helper(curr_k_index, j, pre_next_frontier_cols)

            try:
                _realizing_cols = self._get_realization_templates_dfs_helper(
                    curr_k_index + 1,
                    next_cols,
                    next_frontier_cols
                )

                if len(_realizing_cols) > 0:
                    non_realizable_branch = False
                    realizing_cols.extend(_realizing_cols)

                else:
                    return []

            except Not_Realizable:
                pass

        if non_realizable_branch:
            raise Not_Realizable
        else:
            return realizing_cols

class Realizing_Matrix:
    """This class represents a matrix that realizes a given set.

    This class memoizes its instances, as long as `Realizing_Matrix.remember_instances is True` (the default value is
    `True`). Memoized instances can be accessed via the static method `get_instance`.
    """

    instances = {}

    remember_instances = True

    def __init__(self):
        """This constructor method has no arguments because all attributes of the matrix are set via the methods
        `set_cols` or `set_array`.

        The user will need to explicitly call this constructor method, as `Realizing_Matrix.get_instance` takes a
        `Realizing_Matrix` as input.
        """

        self.array = None
        self.cols = None
        self.K = None
        self.m = None
        self._swap_neighbors = None
        self._swap_class = None
        self._swap_compatibility_matrix = None

    @classmethod
    def get_instance(cls, mat):
        """Get a memoized instance.

        In order to call this method, do something like the following:

            array = np.array( ... ) # replace the ... with your matrix
            max_k = 3 # or whatever you want max_k to be
            mat = Realizing_Matrix()
            mat.set_array(array, max_k)
            mat = Realizing_Matrix.get_instance(mat)

        :param mat: type `Realizing_Matrix`.
        :return: type `Realizing_Matrix`, the memoized instance.
        """

        if not cls.remember_instances:
            return mat

        else:
            if mat in cls.instances.keys():
                return cls.instances[mat]
            else:
                cls.instances[mat] = mat
                return mat

    def set_array(self, array, max_k, skip_not_realizing_matrix_check = False):
        """Set this matrix given a `numpy.ndarray`.

        :param array: (type `numpy.ndarray` with `dtype = int`) The permutation matrix.
        :param max_k: Positive `int`. The is the maximum value of the `Realizable_Set` that this
        matrix realizes.
        :param skip_not_realizing_matrix_check: (type `bool`, default `False`) If `True`, then skip value checks for
        `array`, which may speed-up your code.
        :raises TypeError: If `array` is not a `numpy.ndarray`, or `array` does not have `dtype=int`, or `max_k` is not
        of type `int`.
        :raises Already_Instantiated_Error: If this method or `set_cols` has previously been called on this instance.
        :raises Not_Realizing_Matrix: If `array` does not represent a permutation matrix.
        :raises ValueError: If `max_k` is too small or too large.
        """

        self._check_already_init_raise()

        if (
            not isinstance(array, np.ndarray) or
            not array.dtype == int or
            not isinstance(max_k, (int, np.int32, np.int64))
        ):
            raise TypeError

        m = array.shape[0]

        if (not skip_not_realizing_matrix_check and (
            array.shape == (0,0) or
            array.shape[0] != array.shape[1] or
            np.any(np.sum(array, axis = 0) != 1) or
            len(np.where(array == 1)[0]) != m or
            len(np.where(array == 0)[0]) != m**2 - m
        )):
            raise Not_Realizing_Matrix

        self.array = array
        self.m = m
        self.cols = []
        for i in range(self.m):
            j = np.nonzero(self.array[i, :])[0][0]
            self.cols.append(j)
        self.cols = tuple(self.cols)
        self._set_K(max_k)

    def set_cols(self, cols, max_k, skip_not_realizing_matrix_check = False):
        """Set this matrix given a `list` of column indices. The number `cols[i]` is the column of the `i`-th row
        that is 1; all other columns of the `i`-th row are 0.

        :param cols: A `list` or `numpy.ndarray` of non-negative `int`s.
        :param max_k: Positive `int`.
        :param skip_not_realizing_matrix_check: (type `bool`, default `False`) If `True`, then skip value checks for
        `cols`, which may speed-up your code, but may result in bizarre errors.
        :raises TypeError: If `cols` is not a `list` or a `numpy.ndarray` with `dtype=int`, or if `max_k` is not an
        `int`
        :raises Not_Realizing_Matrix: If `cols` does not represent a permutation matrix.
        :raises Already_Instantiated_Error: If this method or `set_array` has previously been called on this instance.
        :raises ValueError: If `max_k` is too small or too large.
        """

        self._check_already_init_raise()

        if (
            not(isinstance(cols, list) or (isinstance(cols, np.ndarray) and cols.dtype == np.int))
            or not isinstance(max_k, int)
        ):
            raise TypeError

        if (not skip_not_realizing_matrix_check and (
            (isinstance(cols, np.ndarray) and (
                cols.ndim != 1 or
                cols.shape == (0,) or
                np.any(cols < 0) or
                np.any(cols >= len(cols)) or
                len(np.unique(cols)) != len(cols)
            )) or
            (isinstance(cols, list) and (
                len(cols) == 0 or
                not all(isinstance(c, int) for c in cols) or
                any(c < 0 for c in cols) or
                any(c >= len(cols) for c in cols) or
                len(set(cols)) != len(cols)
            ))
        )):
            raise Not_Realizing_Matrix

        self.cols = tuple(cols)
        self.m = len(self.cols)
        self.array = np.zeros((self.m, self.m), dtype=int)
        for i in range(self.m):
            self.array[i, self.cols[i]] = 1
        self._set_K(max_k)

    def _check_init_raise(self, method_name):
        if self.array is None:
            raise Not_Instantiated_Error(method_name)

    def _check_already_init_raise(self):
        if self.array is not None:
            raise Already_Instantiated_Error

    def get_col(self, i):
        return self.cols[i]

    def get_row(self, j):
        return self.cols.index(j)

    def __len__(self):
        self._check_init_raise("len")
        return self.array.shape[0]

    def __hash__(self):
        self._check_init_raise("hash")
        return hash(self.cols) + hash(self.K)

    def __eq__(self, other):
        self._check_init_raise("eq")
        return self.cols == other.cols and self.K == other.K

    def __repr__(self):
        self._check_init_raise("repr")
        return repr(self.array)[:-1] + ", dtype=int)"

    def __str__(self):
        self._check_init_raise("str")
        return "".join(x for x in repr(self) if x in ["0", "1", "\n"])

    def _symmetric_swap(self, i1, i2):
        new_A = np.copy(self.array)
        new_A[i1, :] = self.array[i2, :]
        new_A[i2, :] = self.array[i1, :]
        old_col = np.copy(new_A[:, i1])
        new_A[:, i1] = new_A[:, i2]
        new_A[:, i2] = old_col
        ret = Realizing_Matrix()
        ret.set_array(new_A,self.m-1)
        ret = Realizing_Matrix.get_instance(ret)
        return ret

    def _set_K(self, max_k):

        self._check_init_raise("_set_K")

        if max_k <= 0:
            raise ValueError("`max_k` must be positive.")

        if max_k >= 2*len(self):
            raise ValueError("`max_k` can be at most 2p-1, where p is the size of this matrix")

        upper = np.triu(self.array, len(self) - max_k)
        rows, cols = np.nonzero(upper)
        K = len(self) + rows - cols

        if len(np.unique(K)) != len(K):
            raise Not_Realizing_Matrix

        self.K = Realizable_Set.get_instance(K)

    def delete_one(self, i):
        """Delete the i-th row and j-th column of this matrix, where j is the index of the column of the 1 in the i-th
        row. Does not change this instance.

        :param i: (non-negative `int`) The row to delete.
        :raises Not_Realizing_Matrix: If deleting the i-th row and j-th column gives more than one 1 in any k-diagonal,
        where `k <= self.K.get_max()`.
        :raises TypeError: If `i` is not an `int`.
        :raises ValueError: If `i` is negative.
        :return:
        """

        self._check_init_raise("delete_one")

        if not isinstance(i, int):
            raise TypeError("`i` must be an `int`.")

        if i < 0:
            raise ValueError("`i` must be non-negative.")

        j = self.cols[i]
        ret = np.delete(self.array, i, axis = 0)
        ret = np.delete(ret, j, axis = 1)
        mat = Realizing_Matrix()
        mat.set_array(ret, min(self.K.get_max(), 2*len(self) - 3), False)
        mat = Realizing_Matrix.get_instance(mat)
        return mat

    def get_dims(self):
        upper = np.triu(self.array, 1)
        rows, cols = np.nonzero(upper)
        horz_dim = len(self) - np.min(cols)
        vert_dim = np.max(rows) + 1
        return horz_dim, vert_dim

    def get_swap_neighbors(self):

        self._check_init_raise("get_swap_neighbors")

        if self._swap_neighbors is not None:
            return self._swap_neighbors

        else:
            self._swap_neighbors = []

            for i1, i2 in combinations(range(self.m),2):
                try:
                    self._swap_neighbors.append(self._symmetric_swap(i1,i2))
                except Not_Realizing_Matrix:
                    pass

            return self._swap_neighbors

    def get_swap_class(self):

        self._check_init_raise("get_swap_class")

        if self._swap_class is not None:
            return self._swap_class

        else:
            frontier = {self}
            self._swap_class = {self}
            while len(frontier) > 0:
                new_frontier = set()
                for B in frontier:
                    for neigh in B.get_swap_neighbors():
                        if neigh not in self._swap_class:
                            self._swap_class.add(neigh)
                            new_frontier.add(neigh)
                frontier = new_frontier

            return self._swap_class

    def get_swap_compatibility_matrix(self):

        self._check_init_raise("get_swap_compatibility_matrix")

        if self._swap_compatibility_matrix is not None:
            return self._swap_compatibility_matrix

        else:
            ret = np.zeros((self.m,self.m), dtype=int)
            for i1, i2 in combinations(range(self.m),2):
                try:
                    self._symmetric_swap(i1,i2).cols
                except Not_Realizing_Matrix:
                    continue
                ret[i1, i2] = 1
                ret[i2, i1] = 1


            for i in range(self.m):
                ret[i,i] = 1

            self._swap_compatibility_matrix = ret
            return self._swap_compatibility_matrix

class Realizing_Matrix_Template (Realizing_Matrix):

    instances = {}

    remember_instances = True

    def set_array(self, array, max_k, skip_not_realizing_matrix_check = False):
        """Set this matrix template given a `numpy.ndarray`.

        :param array: (type `numpy.ndarray` with `dtype = int`) 0-1 valued square, strictly upper triangular matrix with
        at most one 1 per column and at most one 1 per row.
        :param max_k: Positive `int`. The is the maximum value of the `Realizable_Set` that this
        matrix template realizes.
        :param skip_not_realizing_matrix_check: (type `bool`, default `False`) If `True`, then skip value checks for
        `array`, which may speed-up your code.
        :raises TypeError: If `array` is not a `numpy.ndarray`, or `array` does not have `dtype=int`, or `max_k` is not
        of type `int`.
        :raises Already_Instantiated_Error: If this method or `set_cols` has previously been called on this instance.
        :raises Not_Realizing_Matrix: If `array` does not represent a realizing matrix template
        :raises ValueError: If `max_k` is too small or too large.
        """

        self._check_already_init_raise()

        if not isinstance(array, np.ndarray) or not array.dtype == int or not isinstance(max_k, int):
            raise TypeError

        m = array.shape[0]

        if (not skip_not_realizing_matrix_check and (
            array.shape == (0,0) or
            array.shape[0] != array.shape[1] or
            np.any((np.sum(array, axis = 0) != 1) * (np.sum(array, axis = 0)) != 0)
        )):
            raise Not_Realizing_Matrix_Template

        elif not skip_not_realizing_matrix_check:
            ones = np.where(array == 1)
            zeroes = np.where(array == 0)
            if np.any(ones[0] >= ones[1]) or len(ones[0]) + len(zeroes[0]) < m ** 2:
                raise Not_Realizing_Matrix_Template

        self.array = array
        self.m = m
        self.cols = [None] * self.m
        for i in range(self.m):
            j = np.nonzero(self.array[i, :])[0]
            if len(j) > 0:
                self.cols[i] = j[0]
        self.cols = tuple(self.cols)
        self._set_K(max_k)

    def set_cols(self, cols, max_k, skip_not_realizing_matrix_check = False):
        """Set this matrix given a `list` of column indices. The number `cols[i]` is the column of the `i`-th row
        that is 1; all other columns of the `i`-th row are 0. If `cols[i] is None`, then all entries of the `i`-th row
        are 0.

        :param cols: A `list` of non-negative `int`s and `None`s
        :param max_k: Positive `int`.
        :param skip_not_realizing_matrix_check: (type `bool`, default `False`) If `True`, then skip value checks for
        `cols`, which may speed-up your code, but may result in bizarre errors.
        :raises TypeError: If `cols` is not a `list`, or if `max_k` is not an `int`
        :raises Not_Realizing_Matrix: If `cols` does not represent a permutation matrix.
        :raises Already_Instantiated_Error: If this method or `set_array` has previously been called on this instance.
        :raises ValueError: If `max_k` is too small or too large.
        """

        self._check_already_init_raise()

        if (
            not isinstance(cols, list) or not isinstance(max_k, int)
        ):
            raise TypeError

        if not skip_not_realizing_matrix_check and (
            len(cols) == 0 or
            not all(isinstance(c, int) or c is None for c in cols) or
            any(not (0 <= c < len(cols)) for c in cols if c is not None) or
            len(set(cols)) != len(cols)
        ):
            raise Not_Realizing_Matrix

        self.cols = tuple(cols)
        self.m = len(self.cols)
        self.array = np.zeros((self.m, self.m), dtype=int)
        for i, c in enumerate(self.cols):
            if c is not None:
                self.array[i, c] = 1
        self._set_K(max_k)

    def iterate_realizations(self):

        all_possible_cols = [
            [j for j in range(i + 1) if j not in self.cols]
            for i in range(len(self))
            if self.cols[i] is None
        ]

        ranges = [
            range(len(cols) - i)
            for i, cols in enumerate(all_possible_cols)
        ]

        last_possible_cols = all_possible_cols[-1]

        for comb in product(*ranges):

            realization_cols = []
            comb_index = 0
            empty = [1] * len(last_possible_cols)

            for i, j in enumerate(self.cols):

                if j is None:
                    for k, acc in enumerate(accumulate(empty)):
                        if acc > comb[comb_index]:
                            break
                        if k >= len(all_possible_cols[comb_index]):
                            raise RuntimeError
                    realization_cols.append(last_possible_cols[k])
                    comb_index += 1
                    empty[k] = 0

                else:
                    realization_cols.append(j)

            mat = Realizing_Matrix()
            mat.set_cols(realization_cols, self.K.get_max(), True)
            mat = Realizing_Matrix.get_instance(mat)
            yield mat

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

def get_maximal_symmetric_swap_distance_matrix(m, cap = 20):

    if m % 2 == 0:
        raise NotImplementedError

    mats = Realizing_Matrix.get_maximal_symmetric_realizing_matrices(m)
    num_mats = len(mats)

    all_neighs = set()
    for mat in mats:
        neighs = mat.get_swap_neighbors()
        for neigh in neighs:
            all_neighs.add((mat, neigh))

    adj_mat = np.zeros((num_mats, num_mats), dtype = int)

    for (mat1, mat2) in all_neighs:
        i1 = mats.index(mat1)
        i2 = mats.index(mat2)
        adj_mat[i1,i2] = 1

    dist_mat = np.copy(adj_mat)
    for i in range(num_mats):
        dist_mat[i,i] = -1

    adj_mat_pow = adj_mat

    d = 1
    while np.any(dist_mat == 0) and d < cap:
        d += 1
        adj_mat_pow = np.matmul(adj_mat_pow, adj_mat)
        for i1, i2 in product(range(num_mats), repeat=2):
            if i1 != i2 and adj_mat_pow[i1,i2] != 0 and dist_mat[i1,i2] == 0:
                dist_mat[i1,i2] = d

    dist_mat[np.nonzero(dist_mat == 0)] = -1
    for i in range(num_mats):
        dist_mat[i,i] = 0

    return dist_mat

def get_symmetric_realizations_(m):

    if m % 2 == 0:
        raise NotImplementedError

    ret = []
    for K in combinations(range(1, m), (m - 1) // 2):
        K = Realizable_Set.get_instance(K)
        if K.check_realizable() > 0:
            ret.extend(K.get_symmetric_realizations(m))
    return ret

def calc_blocking_index(set_length, set_max, check_length, check_max):
    blocking_set_by_realizable_set = {}
    for K in combinations(range(1,set_max+1), set_length):
        K = Realizable_Set.get_instance(K)
        for n in range(check_length + 1):
            for Kp in combinations(range(K[-1]+1, check_max + 1), n):
                Kp = Realizable_Set.get_instance(K.K + Kp)
                if K.is_blocking():
                    blocking_set_by_realizable_set[K] = Kp

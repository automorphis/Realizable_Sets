import copy
import time
from math import floor, ceil
import matplotlib.pyplot as plt

import numpy as np
from itertools import chain, combinations, product

class Not_Realizable (RuntimeError): pass

class Not_Realizing (RuntimeError): pass

class Realizable_Set:

    instances = {}

    remember_instances = True

    @staticmethod
    def get_instance(K):

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
        self.K = tuple(sorted(list(K)))
        self._is_realizable_nece_cond = None
        self._blocks = None
        self._realizations = None
        self._symm_realizing_mats = None
        self._realizable_suff_cond = None
        self._is_realizable = None

    def __len__(self):
        return len(self.K)

    def __iter__(self):
        return iter(self.K)

    def __getitem__(self, index):
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
        return self.K.__format__(format_spec)

    def is_realizable_necessary_condition(self):

        if self._is_realizable_nece_cond is not None:
            return self._is_realizable_nece_cond

        else:

            for i, k in enumerate(self):
                Kp = Realizable_Set.get_instance(self.K[:i+1])
                if (
                    (3/2)*i + 1 > k or (
                        i > 0 and
                        Realizable_Set.get_instance(self.K[:i]).is_blocking() and
                        not self.blocking_condition(i)
                    )
                ):
                    self._is_realizable_nece_cond = False
                    break

            else:
                self._is_realizable_nece_cond = True

            return self._is_realizable_nece_cond

    # def is_realizable_sufficient_condition(self):
    #
    #     if self._realizable_suff_cond is not None:
    #         return self._realizable_suff_cond
    #
    #     else:
    #         for i in range(1, len(self.K)):
    #             Kp = Realizable_Set.get_instance(self.K[:i])
    #             if (
    #                 not Kp.is_realizable() or
    #                 (Kp.is_blocking() and not self.blocking_condition(i)) or
    #                 2*i - 1 >= self.K[i]
    #             ):
    #                 self._realizable_suff_cond = False
    #                 break
    #
    #         else:
    #             self._realizable_suff_cond = True
    #
    #         return self._realizable_suff_cond

    def is_realizable(self):
        if self._is_realizable is not None:
            return self._is_realizable

        else:

            if len(self) == 1:
                self._is_realizable = True
                return self._is_realizable

            if Realizable_Set.remember_instances:
                realizable = True
                for i in range(1, len(self)):
                    Kp = Realizable_Set.get_instance(self.K[:i+1])

                    if not realizable:
                        Kp._is_realizable = False

                    elif Kp._is_realizable is None:
                            Ki = self[i]
                            Kp._is_realizable = Kp._is_realizable_dfs_helper(
                                0,
                                [None] * Ki,
                                list(range(Ki - self[0], Ki))
                            )
                            if not Kp._is_realizable:
                                realizable = False

                    elif not Kp._is_realizable:
                        realizable = False

                return self._is_realizable

            else:
                maxK = self.K[-1]
                self._is_realizable = self._is_realizable_dfs_helper(
                    0,
                    [None] * maxK,
                    list(range(maxK - self[0], maxK))
                )
                return self._is_realizable

    def get_realization_necessary_condition_vector(self):
        return [k - ((3 / 2) * i + 1) for i,k in enumerate(self)]

    def get_realization_necessary_condition_direction_vector(self):
        vect = self.get_realization_necessary_condition_vector()
        return ["u" if v0 <= v1 else "d" for v0,v1 in zip(vect[:-1],vect[1:])]

    def get_realizations(self):

        n = self.K[-1]

        if self._realizations is not None:
            return self._realizations

        else:
            try:
                self._realizations = self._get_realizations_dfs_helper(
                    0,
                    [None]*n,
                    list(range(n-self[0], n))
                )
            except Not_Realizable:
                self._realizations = []

            return self._realizations

    def is_blocking(self):
        return len(self.get_blocks()) > 0

    def get_blocks(self):

        if self._blocks is not None:
            return self._blocks

        else:

            if not self.is_realizable():
                self._blocks = []
                return self._blocks

            elif len(self) < (1 + max(self)) / 2.:
                self._blocks = []

            else:
                m = len(self)
                self._blocks = []
                try:
                    realizing_cols = self._get_blocks_dfs_helper(
                        0,
                        [None]*m,
                        list( range(m-self[0], m) )
                    )

                    for cols in realizing_cols:
                        A = Realizing_Matrix()
                        A.set_cols(cols,self[-1])
                        A = Realizing_Matrix.get_instance(A)
                        self._blocks.append(A)

                except Not_Realizable:
                    pass

            return self._blocks

    def blocking_condition(self, m):
        max_Kp = self[m-1]
        return all(i not in self for i in range(max_Kp + 1, 2 * m + 1))

    def get_symmetric_realizing_matrices(self, m):

        if self._symm_realizing_mats is not None:
            return self._symm_realizing_mats

        else:

            try:
                realizing_cols = self._get_symmetric_realizing_matrices_dfs_helper(
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
                A.set_cols(cols,m-1)
                A = Realizing_Matrix.get_instance(A)
                self._symm_realizing_mats.append(A)

            return self._symm_realizing_mats

    def _is_realizable_dfs_helper(self, curr_k_index, cols, frontier_cols):

        m = self.K[-1]

        if len(frontier_cols) == 0:
            return sum(j is not None for j in cols) == len(self)

        curr_k = self[curr_k_index]
        next_k_index = curr_k_index + 1

        if next_k_index < len(self):
            next_k = self[next_k_index]
            try:
                pre_next_frontier_cols = [
                    j
                    for j in range( max(0, m-next_k),  min(m, 2*m - next_k) )
                    if j not in cols and cols[j + next_k - m] is None
                ]
            except IndexError:
                pass

        for j in frontier_cols:

            i = j - m + curr_k
            next_cols = copy.copy(cols)
            next_cols[i] = j

            if next_k_index < len(self):
                next_frontier_cols = copy.copy(pre_next_frontier_cols)
                for jp in [j, m-next_k+i]:
                    try:
                        next_frontier_cols.remove(jp)
                    except ValueError:
                        pass
            else:
                next_frontier_cols = []

            if self._is_realizable_dfs_helper(
                next_k_index,
                next_cols,
                next_frontier_cols
            ):
                return True

        return False

    def _get_blocks_dfs_helper(self, curr_k_index, cols, frontier_cols):

        if len(frontier_cols) == 0:
            if all(j is not None for j in cols):
                return [cols]
            else:
                raise Not_Realizable

        curr_k = self[curr_k_index]
        next_k_index = curr_k_index + 1

        m = len(self)

        if next_k_index < m:
            next_k = self[next_k_index]
            pre_next_frontier_cols = [
                j
                for j in range( max(0, m-next_k),  min(m, 2*m - next_k))
                if j not in cols and cols[j + next_k - m] is None
            ]

        non_realizable_branch = True
        realizing_cols = []

        for j in frontier_cols:

            i = j - m + curr_k
            next_cols = copy.copy(cols)
            next_cols[i] = j

            if next_k_index < m:
                next_frontier_cols = copy.copy(pre_next_frontier_cols)
                for jp in [j, m-next_k+i]:
                    try:
                        next_frontier_cols.remove(jp)
                    except ValueError:
                        pass
            else:
                next_frontier_cols = []

            try:
                _realizing_cols = self._get_blocks_dfs_helper(
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

    def _get_symmetric_realizing_matrices_dfs_helper(self, m, curr_k_index, cols, frontier_cols):

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
                _realizing_cols = self._get_symmetric_realizing_matrices_dfs_helper(
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

    def _get_realizations_dfs_helper(self, curr_k_index, cols, frontier_cols):

        m = self.K[-1]

        if len(frontier_cols) == 0:
            if sum(j is not None for j in cols) == len(self):
                return [cols]
            else:
                raise Not_Realizable

        curr_k = self[curr_k_index]
        next_k_index = curr_k_index + 1

        if next_k_index < len(self):
            next_k = self[next_k_index]
            pre_next_frontier_cols = [
                j
                for j in range( max(0, m-next_k),  min(m, 2*m - next_k) )
                if j not in cols and cols[j + next_k - m] is None
            ]

        non_realizable_branch = True
        realizing_cols = []

        for j in frontier_cols:

            i = j - m + curr_k
            next_cols = copy.copy(cols)
            next_cols[i] = j

            if next_k_index < len(self):
                next_frontier_cols = copy.copy(pre_next_frontier_cols)
                for jp in [j, m-next_k+i]:
                    try:
                        next_frontier_cols.remove(jp)
                    except ValueError:
                        pass
            else:
                next_frontier_cols = []

            try:
                _realizing_cols = self._get_realizations_dfs_helper(
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

class Realizing_Matrix:

    instances = {}

    remember_instances = True

    def __init__(self):

        self.array = None
        self.cols = None
        self.K = None
        self.m = None
        self._swap_neighbors = None
        self._swap_class = None
        self._swap_compatibility_matrix = None

    @staticmethod
    def get_instance(mat):

        if not Realizing_Matrix.remember_instances:
            return mat

        else:
            if mat in Realizing_Matrix.instances.keys():
                return Realizing_Matrix.instances[mat]
            else:
                Realizing_Matrix.instances[mat] = mat
                return mat

    @staticmethod
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

    @staticmethod
    def get_maximal_symmetric_realizing_matrices(m):

        if m % 2 == 0:
            raise NotImplementedError

        ret = []
        for K in combinations(range(1, m), (m - 1) // 2):
            K = Realizable_Set.get_instance(K)
            if K.check_realizable() > 0:
                ret.extend(K.get_symmetric_realizing_matrices(m))
        return ret

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

    def set_array(self, array, max_k):

        if self.array is not None:
            raise ValueError

        self.array = array
        self.m = len(array)
        self.cols = []
        for i in range(self.m):
            j = np.nonzero(self.array[i, :])[0][0]
            self.cols.append(j)
        self.cols = tuple(self.cols)
        self._set_K(max_k)

    def set_cols(self, cols, max_k):

        if self.cols is not None:
            raise ValueError

        self.cols = tuple(cols)
        self.m = len(self.cols)
        self.array = np.zeros((self.m, self.m), dtype=int)
        for i in range(self.m):
            self.array[i, self.cols[i]] = 1
        self._set_K(max_k)

    def _set_K(self, max_k):
        K = []
        for i,j in enumerate(self.cols):
            k = self.m + i - j
            if k <= max_k:
                if k in K:
                    self.array = None
                    self.cols = None
                    self.K = None
                    raise Not_Realizing
                else:
                    K.append(k)
        self.K = Realizable_Set.get_instance(K)

    def get_swap_neighbors(self):

        if self.array is None:
            raise ValueError

        if self._swap_neighbors is not None:
            return self._swap_neighbors

        else:
            self._swap_neighbors = []

            for i1, i2 in combinations(range(self.m),2):
                try:
                    self._swap_neighbors.append(self._symmetric_swap(i1,i2))
                except Not_Realizing:
                    pass

            return self._swap_neighbors

    def get_swap_class(self):

        if self.array is None:
            raise ValueError

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

        if self.array is None:
            raise ValueError

        if self._swap_compatibility_matrix is not None:
            return self._swap_compatibility_matrix

        else:
            ret = np.zeros((self.m,self.m), dtype=int)
            for i1, i2 in combinations(range(self.m),2):
                try:
                    self._symmetric_swap(i1,i2).cols
                except Not_Realizing:
                    continue
                ret[i1, i2] = 1
                ret[i2, i1] = 1


            for i in range(self.m):
                ret[i,i] = 1

            self._swap_compatibility_matrix = ret
            return self._swap_compatibility_matrix

    def __hash__(self):
        if self.cols is not None:
            return hash(self.cols) + hash(self.K)
        else:
            raise ValueError

    def __eq__(self, other):

        if self.array is None or other.array is None:
            raise ValueError

        return self.cols == other.cols and self.K == other.K

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s) + 1))

# def get_As(m, K):
#
#     K = sorted([m - k for k in K])
#
#
# def get_block(K):
#
#     K = sorted(K)
#
#     if len(K) < (1 + max(K))/2.:
#         return False, []
#
#     else:
#         return _get_block_dfs_helper(K,[],[],list(range(1, 1+min(K))))
#
# def block_condition(K, block_index):
#     max_Kp = K[block_index-1]
#     return all(i not in K for i in range(max_Kp + 1, 2 * block_index + 1))

# def _get_semiblock_dfs_helper(K, rows, cols, frontier_rows):
#
#     curr_k_index = len(rows)
#
#     if curr_k_index >= len(K):
#         max_row, max_col = max(rows), max(cols)
#         if (
#             (max_row == len(K) and max_col == len(K)) or
#             (max_row == len(K) + 1 and max_col + 1 == len(K))
#         ):
#             return True, [rows]
#         else:
#             return False, []
#
#     elif len(frontier_rows) == 0:
#         raise Not_Realizable
#
#     curr_k = K[curr_k_index]
#     next_k_index = curr_k_index + 1
#
#     if next_k_index < len(K):
#         next_k = K[next_k_index]
#         pre_next_frontier_rows = [
#             i
#             for i in range(1,next_k+1)
#             if i not in rows and next_k-i not in cols
#         ]
#
#     realizing_rows = []
#     non_realizable_branch = True
#
#     for i in frontier_rows:
#
#         j = curr_k - i
#         next_rows = rows + [i]
#         next_cols = cols + [j]
#
#         if next_k_index < len(K):
#             next_frontier_rows = copy.copy(pre_next_frontier_rows)
#             try:
#                 next_frontier_rows.remove(i)
#             except ValueError: pass
#             try:
#                 next_frontier_rows.remove(next_k - j)
#             except ValueError: pass
#         else:
#             next_frontier_rows = []
#
#         try:
#             _is_semiblocking, _realizing_rows = _get_semiblock_dfs_helper(
#                 K,
#                 next_rows,
#                 next_cols,
#                 next_frontier_rows
#             )
#
#             non_realizable_branch = False
#             if _is_semiblocking:
#                 realizing_rows.extend(_realizing_rows)
#
#         except Not_Realizable:
#             pass
#
#     if non_realizable_branch:
#         raise Not_Realizable
#     else:
#         return len(realizing_rows) > 0, realizing_rows
#
# def get_semiblock(K):
#     return _get_semiblock_dfs_helper(K,[],[],list(range(1, 1+min(K))))
#
# def get_block_or_semiblock(K):
#
#     K = sorted(K)
#
#     if not is_realizable(K):
#         raise Not_Realizable
#
#     for i in range(1,len(K)+1):
#
#         Kp = K[:i]
#         max_Kp = K[i-1]
#
#         _is_blocking, blocks = get_block(Kp)
#
#         if _is_blocking:
#             return "block", blocks
#
#         _is_semiblocking, semiblocks = get_semiblock(Kp)
#
#         if _is_semiblocking and all(i not in K for i in range(max_Kp+1,2*i-1)):
#             return "semiblock", semiblocks
#
#     else:
#         return None, []

# def is_realizable(K):
#     K = sorted(K)
#     for i in range(1,len(K)):
#         Kp = K[:i]
#         if get_block(Kp)[0] and not block_condition(K, i):
#             return False
#     else:
#         return True
#
# def get_impossible(m, max_subset_size = None, excluded = None):
#     _excluded = []
#     if excluded is not None:
#         for sublist in excluded:
#             _excluded.append(tuple(x+m for x in sublist))
#     excluded = _excluded
#     if max_subset_size is None:
#         max_subset_size = (m-1)//2
#     imp = []
#     for K in powerset(range(1, m)):
#         if 0 < len(K) <= max_subset_size:
#             if len(get_As(m,K)) == 0:
#                 for sublist in excluded:
#                     if K[-len(sublist):] == sublist:
#                         break
#                 else:
#                     imp.append(tuple(m-x for x in K))
#     return imp

# def get_K(A):
#     K = set()
#     m = len(A)
#     for i in range(m):
#         j = np.nonzero(A[i,:])[0][0]
#         k = m+i-j
#         K.add(k)
#     return K

# def get_swap_neighbors(A):
#     m = len(A)
#     ret = []
#     for i1, i2 in combinations(range(m),2):
#         new_A = np.copy(A)
#         new_A[i1,:] = A[i2,:]
#         new_A[i2,:] = A[i1,:]
#         old_col = np.copy(new_A[:,i1])
#         new_A[:,i1] = new_A[:,i2]
#         new_A[:,i2] = old_col
#         K = set()
#         for i in range(m):
#             j = np.nonzero(new_A[i,:])[0][0]
#             k = m+i-j
#             if k != m and k in K:
#                 break
#             elif k != m:
#                 K.add(k)
#         else:
#             ret.append(new_A)
#
#     return ret

# def get_swap_class(A):
#     frontier = set(make_tuple(A))
#     cls = set(make_tuple(A))
#     while len(frontier) > 0:
#         new_frontier = set()
#         for B in frontier:
#             for neigh in get_swap_neighbors(B):
#                 if neigh not in cls:
#                     cls.add(neigh)
#                     new_frontier.add(neigh)
#         frontier = new_frontier
#     return cls


# def make_tuple(A):
#     m = len(A)
#     ret = []
#     for i in range(m):
#         j = np.nonzero(A[i,:])[0][0]
#         ret.append(j)
#     return tuple(ret)
#
# def make_array(tup):
#     m = len(tup)
#     A = np.zeros((m,m),dtype=int)
#     for i in range(m):
#         j = tup[i]
#         A[i,j] = 1
#     return A

# print(get_impossible(13, None,
#                      [
#                          (-2,-1),
#                          (-4,-3,-1),
#                          (-4,-3,-2),
#                          (-7, -6, -5, -3, -2),
#                          (-8, -6, -5, -3, -2),
#                          (-10, -9, -6, -5, -3, -2),
#                          (-4, -3, -1),
#                          (-6, -5, -3, -1),
#                          (-8, -7, -5, -3, -1),
#                          (-10, -9, -7, -5, -3, -1),
#                          (-7, -6, -5, -4, -1),
#                          (-8, -6, -5, -4, -1),
#                          (-10, -9, -6, -5, -4, -1),
#                          (-7, -6, -5, -4, -2)
#                      ]))

# print(is_blocking([2,4,5,6,8]))

# print(get_block_or_semiblock([2,4,5]))

# print(is_realizable([4,5,6]))

# print(get_block_or_semiblock([4,5,6]))

# m = 7
#
# nons = []

# for K in combinations(range(1,m), (m-1)//2):
#     try:
#         kind, block = get_block_or_semiblock(K)
#         As = get_As(m,K)
#         if len(As) == 0 and kind is None:
#             nons.append(K)
#         if kind is None:
#             print(K)
#     except Not_Realizable:
#         pass

# m=7
# K = [3,5,6]
#
# print(get_block_or_semiblock(K))
pass

# A = np.array([
#     [0,0,0,0,1,0,0],
#     [0,0,1,0,0,0,0],
#     [0,1,0,0,0,0,0],
#     [0,0,0,0,0,0,1],
#     [1,0,0,0,0,0,0],
#     [0,0,0,0,0,1,0],
#     [0,0,0,1,0,0,0]
# ])

# A = np.array([
#     [0,1,0],
#     [1,0,0],
#     [0,0,1]
# ])

# A = np.array([
#     [0,1,0,0,0],
#     [1,0,0,0,0],
#     [0,0,0,0,1],
#     [0,0,0,1,0],
#     [0,0,1,0,0]
# ])

# A = np.array([
#     [0,0,0,0,1,0,0,0,0,0,0],
#     [0,0,0,1,0,0,0,0,0,0,0],
#     [0,0,1,0,0,0,0,0,0,0,0],
#     [0,1,0,0,0,0,0,0,0,0,0],
#     [1,0,0,0,0,0,0,0,0,0,0],
#     [0,0,0,0,0,0,0,0,0,0,1],
#     [0,0,0,0,0,0,0,0,0,1,0],
#     [0,0,0,0,0,0,0,0,1,0,0],
#     [0,0,0,0,0,0,0,1,0,0,0],
#     [0,0,0,0,0,0,1,0,0,0,0],
#     [0,0,0,0,0,1,0,0,0,0,0]
# ])

# A = np.array([
#     [0,0,0,0,0,0,0,0,0,0,1],
#     [0,0,0,0,0,0,0,0,0,1,0],
#     [0,0,0,0,0,0,0,0,1,0,0],
#     [0,0,0,0,0,0,0,1,0,0,0],
#     [0,0,0,0,0,0,1,0,0,0,0],
#     [0,0,0,0,0,1,0,0,0,0,0],
#     [0,0,0,0,1,0,0,0,0,0,0],
#     [0,0,0,1,0,0,0,0,0,0,0],
#     [0,0,1,0,0,0,0,0,0,0,0],
#     [0,1,0,0,0,0,0,0,0,0,0],
#     [1,0,0,0,0,0,0,0,0,0,0]
# ])

# mat = Realizing_Matrix()
# mat.set_array(A)

# swap_Ks = set(A.K for A in mat.get_swap_neighbors())

# cls = mat.get_swap_class()
# symm = get_number_maximal_symmetric_realizing_matrices(11)

# x = mat.get_swap_compatibility_matrix()
#
pass


# mat.set_cols(mat.cols)
# pass

# K = Realizable_Set([2,3,4,7])
# print(K.check_realizable())

# get_swap_class(A)

# print(np.all(make_array(make_tuple(A)) == A))
#
# neigh_Ks = set()
#
# for neigh in get_swap_neighbors(A):
#     neigh_Ks.add(tuple(sorted(list(get_K(neigh)))[:3]))
#
# print(neigh_Ks)
#
# print(get_block([2,5,6]))

# print(Realizable_Set([2,3]).get_blocks())

# print(Realizable_Set([2,3,4]).check_realizable())

# mat = Realizing_Matrix()
# mat.set_array(A, len(A)-1)
# mat = Realizing_Matrix.get_instance(mat)
#
# cls1 = mat.get_swap_class()
# cls2 = Realizing_Matrix.get_maximal_symmetric_realizing_matrices(len(A))

# for mat in cls2:
#     if cls2.count(mat) > 1 or len(mat.K) < 3:
#         pass

# mat = Realizing_Matrix.get_maximal_symmetric_swap_distance_matrix(11)
# print(np.max(mat))

# print(Realizable_Set([2,4,5,6,7]).is_realizable_sufficient_condition())

pass

# print(Realizable_Set([3,4,5,6,7]).get_realizations())

# mat = Realizing_Matrix()
# mat.set_cols([4,2,5,3,6,0,1],7)
pass

# print(Realizable_Set.get_instance([6,7,8,9,10,11,12,13,14,15,16]).is_blocking())

def calc_blocking_index(set_length, set_max, check_length, check_max):
    blocking_set_by_realizable_set = {}
    for K in combinations(range(1,set_max+1), set_length):
        K = Realizable_Set.get_instance(K)
        for n in range(check_length + 1):
            for Kp in combinations(range(K[-1]+1, check_max + 1), n):
                Kp = Realizable_Set.get_instance(K.K + Kp)
                if K[.is_blocking():
                    blocking_set_by_realizable_set[K] = Kp





# print(Realizable_Set.get_instance([2,4,5,6,7]).is_realizable_necessary_condition())
#
# Realizable_Set.remember_instances = True
# Realizing_Matrix.remember_instances = True
#
# maxK = 10
# maxN = 7
# realX = []
# realY = []
# counterX = []
# counterY = []
# for n in range(1, maxN+1):
#     start = time.time()
#     print(n)
#     for K in combinations(range(1,maxK+1), n):
#         K = Realizable_Set.get_instance(K)
#         vect = K.get_realization_necessary_condition_vector()
#         dirs = K.get_realization_necessary_condition_direction_vector()
#         if K.is_realizable():
#             for i,k in enumerate(K):
#                 realX.append(i+1)
#                 realY.append(k)
#             if any(d == 0 or d == 1/2 for d in vect):
#                 print(f"good {K} {vect} {''.join(dirs)} {sum(vect)/len(K)}")
#         elif K.is_realizable_necessary_condition():
#             print(f"bad  {K} {vect} {''.join(dirs)} {sum(vect)/len(K)}")
#     print(time.time() - start)
#
#     # print(N,1, len(counterexamples1))
#     # print(N,2, len(counterexamples2))
#
# start = time.time()
# plt.figure()
# plt.scatter(realX,realY, s=20, c="r")
# plt.scatter(counterX, counterY, s=10, c="b")
# I = list(range(1,max(realX)+1))
# plt.plot(I, [ceil((3/2)*i - 0.5) for i in I])
# # plt.plot(I, [2*i - 3 for i in I])
# print(time.time() - start)
# plt.show()




# x = K.get_symmetric_realizing_matrices(3)

import time

from realize import Realizable_Set, powerset, Realizing_Matrix_Template, Realizing_Matrix
from itertools import combinations, product


def realizable_sets(subset_sizes, max_elements):
    for subset_size, max_element in product(subset_sizes, max_elements):
        if subset_size <= max_element:
            if subset_size == 1:
                yield Realizable_Set.get_instance((max_element,))
            else:
                for subset in combinations(range(1,max_element), subset_size-1):
                    K = Realizable_Set.get_instance(subset + (max_element,))
                    if K.is_realizable():
                        yield K

Realizing_Matrix_Template.remember_instances = False
Realizable_Set.remember_instances = False
Realizing_Matrix.remember_instances = False

# any blocking realization => blocking conjecture
for subset_size, max_element in product(range(9, 12), range(1, 18)):
    if 2*subset_size > max_element >= subset_size:
        start = time.time()
        num_blocking = 0
        for K in realizable_sets([subset_size], [max_element]):
            _, fat_dim = K.get_skinny_fat_rectangle_dimensions()
            if fat_dim[0] == fat_dim[1] and fat_dim[0] == len(K) and not K.is_blocking():
                print(1, K)
            elif K.is_blocking():
                num_blocking += 1
                if fat_dim[0] != fat_dim[1] or fat_dim[0] != len(K):
                    print(2, K)
        print(subset_size, max_element, num_blocking, time.time() - start)

# K\cap [1, 2i] = K' conjecture
# for subset_size, max_element in product(range(1, 12), range(1, 18)):
#     start = time.time()
#     for K in realizable_sets([subset_size], [max_element]):
#         for Kp in powerset(K.K):
#             if 0 < len(Kp) < len(K):
#                 Kp = Realizable_Set.get_instance(Kp)
#                 i = K.index(Kp.get_max()) + 1
#                 if Kp.is_blocking() and len([j for j in range(1, 2*i+1) if j in K]) != len(Kp):
#                     print(K,Kp)
#     print(subset_size, max_element, time.time() - start)

# delete non-maximal element conjecture
# for subset_size, max_element in product(range(2, 12), range(1, 18)):
#     start = time.time()
#     for K in realizable_sets([subset_size], [max_element]):
#         if K.is_blocking():
#             for Kp in combinations(K.K, subset_size-1):
#                 if max(Kp) == max_element:
#                     if Realizable_Set.get_instance(Kp).is_blocking():
#                         print(Kp, K)
#     print(subset_size, max_element, time.time() - start)

# subtract 1 from blocking set conjecture
for subset_size, max_element in product(range(2, 12), range(1, 18)):
    start = time.time()
    for K in realizable_sets([subset_size], [max_element]):
        if K.is_blocking():
            for i, k in enumerate(K):
                if k-1 not in K:
                    Kp = Realizable_Set.get_instance(K.K[:i] + K.K[i+1:])
                    if Kp.is_realizable():
                        print(K,Kp)
    print(subset_size, max_element, time.time() - start)

# for subset_size, max_element in product(range(1, 12), range(1, 18)):


# print(Realizable_Set.get_instance((3,)).get_skinny_fat_rectangle_dimensions())

# n = 7
# print(all(Realizable_Set.get_instance(range(n,k+1)).is_blocking() for k in range(n,3*n-2)))

# print(len(Realizable_Set.get_instance(range(n,3*n-1)).get_blocks()))
# x=None

# n = 8
# print(Realizable_Set.get_instance(range(n,3*n-1)).is_blocking())
# print(Realizable_Set.get_instance((n-1,) + tuple(range(n+1,3*n-1))).is_realizable())
#
# n = 9
# print(Realizable_Set.get_instance(range(n,3*n-1)).is_blocking())
# print(Realizable_Set.get_instance((n-1,) + tuple(range(n+1,3*n-1))).is_realizable())

# largest_element = 20
# largest_subset_size = 10
#
# ret = []
#
# for subset_size in range(1, largest_subset_size + 1):
#     print(subset_size)
#     for comb in combinations(range(1,largest_element+1), subset_size):
#         rs = Realizable_Set.get_instance(comb)
#         if rs.is_realizable():
#             for rsp in powerset(rs):
#                 if 0 < len(rsp) < len(rs):
#                     rsp = Realizable_Set.get_instance(rsp)
#                     k0 = max(rsp)
#                     i = rsp.K.index(k0) + 1
#                     if rsp.is_blocking() and [j for j in range(1,2*i+1) if j in rs] != list(rsp):
#                         ret.append((rs, rsp))
#                         print(rs, rsp)

# for x in ret:
#     print(x)


# largest = 20
# set_size = 8
#
# ret = []
#
# for max_k in range(1, largest + 1):
#     for comb in combinations(range(1,max_k), set_size-1):
#         comb += (max_k,)
#         st1 = Realizable_Set.get_instance(comb)
#         if st1.is_blocking():
#             for j in range(min(comb)+1,max(comb)):
#                 if j not in st1:
#                     st2 = Realizable_Set.get_instance(comb + (j,))
#                     if st2.is_realizable():
#                         ret.append((st1,j))
#
# for st, num in ret:
#     print(st, num)

# largest = 10
# set_size = 4
#
# ret = []
#
# for max_k in range(1, largest + 1):
#     for comb in combinations(range(1,max_k), set_size-1):
#         comb += (max_k,)
#         st = Realizable_Set.get_instance(comb)
#         if not st.is_realizable():
#             ret.append(st)
#
# for st in ret:
#     print(st)
#
# print(Realizable_Set.get_instance((2,4,5)).is_blocking())
# print(Realizable_Set.get_instance((2,4,5,6)).is_blocking())

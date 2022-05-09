from main import Realizable_Set, powerset
from itertools import combinations

print(Realizable_Set.get_instance((3,)).get_skinny_fat_rectangle_dimensions())

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

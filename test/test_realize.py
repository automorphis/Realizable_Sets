from itertools import product, combinations
from unittest import TestCase

import numpy as np

from realize import Realizable_Set, Realizing_Matrix, Realizing_Matrix_Template

realizable_sets = [
    (1, ),
    (1, 3),
    (1, 4),
    (1, 3, 5),
    (1, 3, 6),
    (2, ),
    (2, 3),
    (2, 3, 4),
    (2, 3, 4, 7),
    (2, 3, 4, 8),
    (2, 4),
    (2, 4, 5),
    (2, 4, 5, 6),
    (2, 4, 5, 6, 8),
    (2, 4, 5, 6, 8, 11),
    (3, ),
    (3, 4),
    (3, 4, 5),
    (3, 4, 5, 6),
    (3, 4, 5, 6, 7),
    (3, 4, 5, 6, 7, 11)
]

non_realizable_sets = [
    (1, 2),
    (1, 2, 3),
    (1, 3, 4),
    (2, 3, 4, 5),
    (2, 3, 4, 6),
    (2, 4, 5, 6, 7),
    (2, 4, 5, 6, 8, 9),
    (2, 4, 5, 6, 8, 10),
    (3, 4, 5, 6, 7, 8),
    (3, 4, 5, 6, 7, 9),
    (3, 4, 5, 6, 7, 10)
]

blocks = [
    ((1,), [

        [[1]]
    ]),

    ((1, 3), [

        [[0, 1],
         [1, 0]]
    ]),

    ((1, 3, 5), [

        [[0, 0, 1],
         [0, 1, 0],
         [1, 0, 0]]
    ]),

    ((2, 3, 4), [

        [[0, 1, 0],
         [1, 0, 0],
         [0, 0, 1]],

        [[1, 0, 0],
         [0, 0, 1],
         [0, 1, 0]]
    ]),

    ((2, 3, 5, 6), [

        [[0, 0, 1, 0],
         [1, 0, 0, 0],
         [0, 0, 0, 1],
         [0, 1, 0, 0]],

        [[0, 1, 0, 0],
         [0, 0, 0, 1],
         [1, 0, 0, 0],
         [0, 0, 1, 0]]
    ]),

    ((2, 4, 5, 6, 8), [

        [[0, 0, 0, 1, 0],
         [0, 0, 1, 0, 0],
         [0, 1, 0, 0, 0],
         [1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1]],

        [[1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 1, 0],
         [0, 0, 1, 0, 0],
         [0, 1, 0, 0, 0]],

        [[0, 0, 0, 1, 0],
         [1, 0, 0, 0, 0],
         [0, 0, 1, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 1, 0, 0, 0]],

        [[0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 1, 0, 0],
         [1, 0, 0, 0, 0],
         [0, 0, 0, 1, 0]]
    ]),

    ((3, 4, 5, 6, 7), [

        [[0, 0, 1, 0, 0],
         [1, 0, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1]],

        [[0, 0, 1, 0, 0],
         [0, 1, 0, 0, 0],
         [1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 1, 0]],

        [[1, 0, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 1, 0, 0]],

        [[1, 0, 0, 0, 0],
         [0, 0, 1, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 1, 0, 0, 0],
         [0, 0, 0, 1, 0]],

        [[0, 1, 0, 0, 0],
         [1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 1, 0],
         [0, 0, 1, 0, 0]],

        [[0, 1, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [1, 0, 0, 0, 0],
         [0, 0, 1, 0, 0],
         [0, 0, 0, 0, 1]]
    ])
]

realizable_non_blocks = [
    (2, ),
    (2, 3),
    (1, 4),
    (1, 5),
    (1, 6),
    (2, 4),
    (2, 5),
    (3, ),
    (3, 4),
    (3, 4, 5),
    (3, 4, 5, 6),
    (4, 5, 6, 7, 8, 9)
]

non_blocking_realization_templates = [

    ((1, 4), [

       [[0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]]
    ]),

    ((1, 5), [

       [[0, 0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]],

       [[0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

    ]),

    ((2 ,), [

        [[0, 1, 0],
         [0, 0, 0],
         [0, 0, 0]],

        [[0, 0, 0],
         [0, 0, 1],
         [0, 0, 0]]
    ]),

    ((2, 3), [

        [[0, 0, 1, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 1],
         [0, 0, 0, 0]],

        [[0, 1, 0, 0],
         [0, 0, 0, 1],
         [0, 0, 0, 0],
         [0, 0, 0, 0]]
    ]),

    ((2, 4), [

        [[0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0]],

        [[0, 0, 0, 1, 0],
         [0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]],

        [[0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]]
    ]),

    ((3, ), [

        [[0, 1, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0]],

        [[0, 0, 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 0]],

        [[0, 0, 0, 0],
         [0, 0, 0, 0],
         [0, 0, 0, 1],
         [0, 0, 0, 0]]
    ]),

    ((3, 4), [

        [[0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]],

        [[0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0]],

        [[0, 1, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0]],

        [[0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0],
         [0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0]]
    ]),

    ((3, 4, 5), [

        [[0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 1, 0, 0],
         [0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

        [[0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

        [[0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

        [[0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

        [[0, 1, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

        [[0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 1, 0, 0],
         [0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0]],
    ]),

    ((3, 4, 5, 6), [

        [[0, 0, 0, 0, 1, 0, 0],
         [0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 0, 0, 1, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],

        [[0, 1, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0]],
    ])

]

non_blocking_skinny_fat_rectangle_dimensions = [
    ((2, ), (1, 2), (1, 2)),
    ((2, 3), (2, 3), (2, 3)),
    ((1, 4), (2, 3), (2, 3)),
    ((1, 5), (2, 4), (3, 3)),
    ((1, 6), (2, 5), (3, 4)),
    ((2, 4), (2, 4), (2, 3)),
    ((2, 5), (2, 5), (3, 3)),
    ((3, ), (1, 3), (2, 2)),
    ((3, 4), (2, 4), (3, 3)),
    ((3, 4, 5), (3, 5), (4, 4)),
    ((3, 4, 5, 6), (4, 6), (5, 5))
]

class Test_Realizable_Set(TestCase):

    def setUp(self):
        Realizable_Set.remember_instances = True
        Realizing_Matrix.remember_instances = True
        Realizing_Matrix_Template.remember_instances = True

    def test_is_realizable(self):

        for _ in range(2):

            for K in realizable_sets:
                self.assertTrue(Realizable_Set.get_instance(K).is_realizable())

            for K in non_realizable_sets:
                self.assertFalse(Realizable_Set.get_instance(K).is_realizable())

            Realizable_Set.remember_instances = False
            Realizing_Matrix.remember_instances = False
            Realizing_Matrix_Template.remember_instances = False

    def test_get_blocks(self):

        self.assertTrue(all(
            all(exp_mat1 != exp_mat2 for exp_mat1, exp_mat2 in combinations(exp_blocks, 2))
            for _, exp_blocks in blocks
        ))

        for _ in range(2):

            for K, exp_blocks in blocks:

                K = Realizable_Set.get_instance(K)
                calc_blocks = K.get_blocks()

                self.assertEqual(len(calc_blocks), len(exp_blocks))

                for exp_block in exp_blocks:

                    exp_mat = Realizing_Matrix()
                    exp_mat.set_array(np.array(exp_block, dtype=int), K.get_max())
                    exp_mat = Realizing_Matrix.get_instance(exp_mat)
                    self.assertTrue(any(exp_mat == mat for mat in calc_blocks))

            for K in realizable_non_blocks:
                K = Realizable_Set.get_instance(K)
                self.assertEqual(K.get_blocks(), [])

            Realizable_Set.remember_instances = False
            Realizing_Matrix.remember_instances = False
            Realizing_Matrix_Template.remember_instances = False

    def test_get_realization_templates(self):

        self.assertTrue(all(
            all(exp_mat1 != exp_mat2 for exp_mat1, exp_mat2 in combinations(exp_real_temps, 2))
            for _, exp_real_temps in non_blocking_realization_templates
        ))

        for _ in range(2):

            for K, exp_blocks in blocks:

                K = Realizable_Set.get_instance(K)

                kp = K.get_max()
                size = kp + 1
                n = len(K)

                calc_real_temp_mats = K.get_realization_templates()

                self.assertEqual(len(calc_real_temp_mats), len(exp_blocks))

                for exp_block in exp_blocks:
                    exp_real_temp = np.pad(exp_block, [(0, size - n), (size - n, 0)])
                    exp_real_temp_mat = Realizing_Matrix_Template()
                    exp_real_temp_mat.set_array(exp_real_temp, kp)
                    exp_real_temp_mat = Realizing_Matrix_Template.get_instance(exp_real_temp_mat)
                    self.assertTrue(any(exp_real_temp_mat == calc_real_temp_mat for calc_real_temp_mat in calc_real_temp_mats))

            for K, exp_real_temps in non_blocking_realization_templates:

                K = Realizable_Set.get_instance(K)

                kp = K.get_max()

                calc_real_temp_mats = K.get_realization_templates()

                self.assertEqual(len(calc_real_temp_mats), len(exp_real_temps))

                for exp_real_temp in exp_real_temps:

                    exp_real_temp = np.array(exp_real_temp)
                    exp_real_temp_mat = Realizing_Matrix_Template()
                    exp_real_temp_mat.set_array(exp_real_temp, kp)
                    exp_real_temp_mat = Realizing_Matrix_Template.get_instance(exp_real_temp_mat)
                    self.assertTrue(any(exp_real_temp_mat == calc_real_temp_mat for calc_real_temp_mat in calc_real_temp_mats))

            Realizable_Set.remember_instances = False
            Realizing_Matrix.remember_instances = False
            Realizing_Matrix_Template.remember_instances = False

    def test_get_skinny_fat_rectangle_dimensions(self):

        for _ in range(2):

            for K, exp_blocks in blocks:

                K = Realizable_Set.get_instance(K)
                n = len(K)
                self.assertEqual(K.get_skinny_fat_rectangle_dimensions(), ((n, n), (n, n)))

            for K, skinny_dims, fat_dims in non_blocking_skinny_fat_rectangle_dimensions:

                self.assertEqual(
                    Realizable_Set.get_instance(K).get_skinny_fat_rectangle_dimensions(),
                    (skinny_dims, fat_dims)
                )

            Realizable_Set.remember_instances = False
            Realizing_Matrix.remember_instances = False
            Realizing_Matrix_Template.remember_instances = False

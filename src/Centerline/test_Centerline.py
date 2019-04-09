#!/usr/bin/env python
import os
import pytest
import numpy as np

import Centerline

class TestCenterLine():
    @pytest.fixture(scope='class')
    def centerline_tester(self):
        xx = np.arange(10)/10.0
        yy = np.arange(10)/10.0
        return Centerline.Centerline(xx, yy)

    def test_distance(self, centerline_tester):
        i, d, x, y, s, n = centerline_tester(0.1, 100)
        assert d == pytest.approx(99.10322901, abs=0.001)

    def test_along_reach(self, centerline_tester):
        i, d, x, y, s, n = centerline_tester(0.1, 100)
        assert s == pytest.approx(69.50859659, abs=0.001)

    def test_cross_reach(self, centerline_tester):
        i, d, x, y, s, n = centerline_tester(0.1, 100)
        assert n == pytest.approx(70.63996744, abs=0.001)

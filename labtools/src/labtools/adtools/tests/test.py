import unittest
import finder

class TestReadMethods(unittest.TestCase):
    """Test the methods for raw read parsing."""

    def setUp(self):
        self.normal_read = "TACTGCGGGCTCTACTTCATCGGCTAGCATTGCTAATTATGCTCATTT"\
        "GGGTCCACAGAACTTTAAGAGAAATTTGCAACAATGTCAGAAGATTGTCTTGAATCCATCTAATATT"\
        "CAATTGAATACTCCACCACAATTTAGATTGGCTTGATAACTAGCTGAGGGCCCGCTACTTCTCCAGG"\
        "CGCGAAACTTATAAGTAAGCGAA"
        self.truncated_read ="TACTGCGGGCTCTACTTCATCGGCTAGCATTGCTAATTATG"\
        "CTCATTTGGGTCCACAGAACTTTAAGAGAAATTTGCAACAATGTCAGAAGATTGTCTTGAATCCATC"\
        "TAATATTCAATTGAATACTCCACCACAATTTAGATTGGCTTGATAACTAGCTG"
        self.random_read ="AAAAAGAGAGAGAGAGAGAGAATGATGCGCGCGCGCTGAGCTAGCTGACTGCT"\
        "AGTAGGCTGATCGATGACGTAGCTAGCTACGACTGCATGACTCAGCATCGCAAAATGCGCGTGGGCT"\
        "GGGGTCGATCGAGTGCATGACTAGCATCAGTCGGCTAGTCGACTACGACTCGATCAGCTAGCAGCAT"\
        "AAGTGGGGCTAGCTAGCGGATCGCATG"

    def test_tile(self):
        """Test the pull_AD function with normal input."""

        expect_AD = "ATTGCTAATTATGCTCATTTGGGTCCACAGAACTTTAAGAGAAATTTGCAACAAT"\
        "GTCAGAAGATTGTCTTGAATCCATCTAATATTCAATTGAATACTCCACCACAATTTAGATTGGCT"
        expect_bc = "CTACTTCTCCA"
        self.assertEqual(finder.pull_AD(self.normal_read)[0], expect_AD)
        self.assertEqual(finder.pull_AD(self.normal_read)[1], expect_bc)
    
    def test_truncated_read(self):
        """Test the pull_AD function with a truncated read."""

        self.assertEqual(finder.pull_AD(self.truncated_read)[0], None)
        self.assertEqual(finder.pull_AD(self.truncated_read)[1], None)

    def test_random_read(self):
        """Test the pull_AD function with an undesigned read."""

        self.assertEqual(finder.pull_AD(self.random_read)[0], None)
        self.assertEqual(finder.pull_AD(self.random_read)[1], None)


class TestExperimentMethods(unittest.TestCase):
    """Tests for ExperimentBuilder Module."""

    pass

class TestSequenceMethods(unittest.TestCase):
    """Tests for Sequence methods."""

    def test_AD_list_equality(self):
        pass

if __name__ == "__main__":

    unittest.main()
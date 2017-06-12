noinst_LIBRARIES = libgtest.a
libgtest_a_SOURCES = gtest/gtest-all.cc gtest/gtest_main.cc gtest/gtest.h

INCLUDES = -I$(top_srcdir)/src/ $(BOOST_CPPFLAGS) $(GFLAGS_CPPFLAGS)
LDADD = ../src/liblampsearch.a libgtest.a $(GFLAGS_LIBS)
TESTS = database_lcm_unittest database_reduction_unittest database_unittest dfs_unittest graph_unittest \
        lamp_graph_unittest lamp_unittest lcm_dfs_lamp_unittest lcm_dfs_lamp_vba_unittest \
        lcm_dfs_unittest lcm_unittest table_vba_unittest utils_unittest variable_bitset_unittest \
        variable_length_itemset_unittest functions_fisher_unittest functions_chi_unittest functions_utest_unittest

check_PROGRAMS = database_lcm_unittest database_reduction_unittest database_unittest dfs_unittest graph_unittest \
                 lamp_graph_unittest lamp_unittest lcm_dfs_lamp_unittest lcm_dfs_lamp_vba_unittest \
                 lcm_dfs_unittest lcm_unittest table_vba_unittest utils_unittest variable_bitset_unittest \
                 variable_length_itemset_unittest functions_fisher_unittest functions_chi_unittest functions_utest_unittest

database_lcm_unittest_SOURCES = database_lcm_unittest.cc

database_reduction_unittest_SOURCES = database_reduction_unittest.cc

database_unittest_SOURCES = database_unittest.cc

dfs_unittest_SOURCES = dfs_unittest.cc

graph_unittest_SOURCES = graph_unittest.cc

lamp_graph_unittest_SOURCES = lamp_graph_unittest.cc

lamp_unittest_SOURCES = lamp_unittest.cc

lcm_dfs_lamp_unittest_SOURCES = lcm_dfs_lamp_unittest.cc

lcm_dfs_lamp_vba_unittest_SOURCES = lcm_dfs_lamp_vba_unittest.cc

lcm_dfs_unittest_SOURCES = lcm_dfs_unittest.cc

lcm_unittest_SOURCES = lcm_unittest.cc

table_vba_unittest_SOURCES = table_vba_unittest.cc

utils_unittest_SOURCES = utils_unittest.cc

variable_bitset_unittest_SOURCES = variable_bitset_unittest.cc

variable_length_itemset_unittest_SOURCES = variable_length_itemset_unittest.cc

functions_fisher_unittest_SOURCES = functions_fisher_unittest.cc

functions_chi_unittest_SOURCES = functions_chi_unittest.cc

functions_utest_unittest_SOURCES = functions_utest_unittest.cc

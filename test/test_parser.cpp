#include "catch.hpp"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "parser.hpp"

SCENARIO( "parser checks for mutually exclusive options", "[parser]" )
{
    GIVEN( "a call to parse_options" )
    {
        std::stringstream out;

        WHEN( "both --init_size and --init_diversity are specified" )
        {
            char *args[] = {"my_program",
                            "--init_size", "25",
                            "--init_diversity", "foo.csv", "bar.csv"};

            THEN( "a logic error is thrown" )
            {
                REQUIRE_THROWS_AS( parser::parse_options(6, args, out), std::logic_error );
            }
        }

        WHEN( "both --load_snapshot and --save_snapshot are specified" )
        {
            char *args[] = {"my_program",
                            "--load_snapshot", "archive.tar",
                            "--save_snapshot"};

            THEN( "a logic error is thrown" )
            {
                REQUIRE_THROWS_AS( parser::parse_options(4, args, out), std::logic_error );
            }
        }
    }
}

SCENARIO( "parser requires param set and run number to be ints", "[parser]" )
{
    GIVEN( "a call to parse_options" )
    {
        std::stringstream out;

        WHEN( "non-integer --param_set is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--param_set", "ten"};
            char *float_arg[] = {"my_program",
                                 "--param_set", "10.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, float_arg, out) );
            }
        }

        WHEN( "non-integer --run_number is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--run_number", "ten"};
            char *float_arg[] = {"my_program",
                                 "--param_set", "10.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, float_arg, out) );
            }
        }
    }
}

SCENARIO( "parser requires max cycles and population size to be ints", "[parser]")
{
    GIVEN( "a call to parse_options" )
    {
        std::stringstream out;

        WHEN( "non-integer --max_cycles is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--max_cycles", "tenthousand"};
            char *float_arg[] = {"my_program",
                                 "--max_cycles", "10000.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, float_arg, out) );
            }
        }

        WHEN( "non-integer --max_size_lim is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--max_size_lim", "tenthousand"};
            char *float_arg[] = {"my_program",
                                 "--max_size_lim", "10000.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, float_arg, out) );
            }
        }
    }
}

SCENARIO( "parser checks that probability values are valid", "[parser]")
{
    GIVEN( "a call to parse_options" )
    {
        std::stringstream out;

        WHEN( "reasonable value for --prob_mut_pos is specified" )
        {
            char *float_arg[] = {"my_program",
                               "--prob_mut_pos", "0.99"};

            THEN( "no exception is thrown" )
            {
                REQUIRE_NOTHROW( parser::parse_options(3, float_arg, out) );
            }
        }

        WHEN( "unreasonable float value for --prob_mut_pos is specified" )
        {
            char *too_large_arg[] = {"my_program", "--prob_mut_pos", "5.0"};
            char *too_small_arg[] = {"my_program", "--prob_mut_pos", "-3.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, too_large_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, too_small_arg, out) );
            }
        }

        WHEN( "non-float --prob_mut_pos is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--prob_mut_pos", "pretty_high_chance"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
            }
        }

        WHEN( "reasonable value for --prob_mut_neg is specified" )
        {
            char *float_arg[] = {"my_program",
                               "--prob_mut_neg", "0.99"};

            THEN( "no exception is thrown" )
            {
                REQUIRE_NOTHROW( parser::parse_options(3, float_arg, out) );
            }
        }
        WHEN( "unreasonable float value for --prob_mut_neg is specified" )
        {
            char *too_large_arg[] = {"my_program", "--prob_mut_neg", "5.0"};
            char *too_small_arg[] = {"my_program", "--prob_mut_neg", "-3.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, too_large_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, too_small_arg, out) );
            }
        }
        WHEN( "non-float --prob_mut_neg is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--prob_mut_neg", "pretty_high_chance"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
            }
        }
        WHEN( "reasonable value for --prob_inc_mut is specified" )
        {
            char *float_arg[] = {"my_program",
                               "--prob_inc_mut", "0.99"};

            THEN( "no exception is thrown" )
            {
                REQUIRE_NOTHROW( parser::parse_options(3, float_arg, out) );
            }
        }
        WHEN( "unreasonable float value for --prob_inc_mut is specified" )
        {
            char *too_large_arg[] = {"my_program", "--prob_inc_mut", "5.0"};
            char *too_small_arg[] = {"my_program", "--prob_inc_mut", "-3.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, too_large_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, too_small_arg, out) );
            }
        }
        WHEN( "non-float --prob_inc_mut is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--prob_inc_mut", "pretty_high_chance"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
            }
        }
        WHEN( "reasonable value for --prob_dec_mut is specified" )
        {
            char *float_arg[] = {"my_program",
                               "--prob_dec_mut", "0.99"};

            THEN( "no exception is thrown" )
            {
                REQUIRE_NOTHROW( parser::parse_options(3, float_arg, out) );
            }
        }
        WHEN( "unreasonable float value for --prob_dec_mut is specified" )
        {
            char *too_large_arg[] = {"my_program", "--prob_dec_mut", "5.0"};
            char *too_small_arg[] = {"my_program", "--prob_dec_mut", "-3.0"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, too_large_arg, out) );
                REQUIRE_THROWS( parser::parse_options(3, too_small_arg, out) );
            }
        }
        WHEN( "non-float --prob_dec_mut is specified" )
        {
            char *str_arg[] = {"my_program",
                               "--prob_dec_mut", "pretty_high_chance"};

            THEN( "an exception is thrown" )
            {
                REQUIRE_THROWS( parser::parse_options(3, str_arg, out) );
            }
        }
    }
}


#include "catch.hpp"
#include <stdexcept>
#include "parser.hpp"

SCENARIO( "parser checks for mutually exclusive options", "[parser]" )
{
    GIVEN( "a call to parser::parse_options" )
    {
        WHEN( "both --init_size and --init_diversity are specified" )
        {
            char *argv[] = {"--init_size", "25",
                            "--init_diversity", "foo.csv", "bar.csv"};

            THEN( "a logic error is thrown" )
            {
                REQUIRE_THROWS_AS( parser::parse_options(5, argv), std::logic_error );
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Keep this file empty, and implement unit tests in separate compilation units!
////////////////////////////////////////////////////////////////////////////////

// Catch2 Documentation: https://github.com/catchorg/Catch2/tree/master/docs

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

#include <spdlog/spdlog.h>

int main(int argc, char* argv[])
{
    Catch::Session session; // There must be exactly one instance

    int log_level = spdlog::level::off;

    // Build a new parser on top of Catch's
    using namespace Catch::clara;
    auto cli = session.cli()
        | Opt(
            [&log_level](int const d) {
                if (d < 0 || d > spdlog::level::off) {
                    return ParserResult::runtimeError(
                        "Log level must be between 0 and 6");
                }
                else {
                    log_level = d;
                    return ParserResult::ok(ParseResultType::Matched);
                }
            },
            "log_level")["-g"]["--logger-level"](
            "logger verbosity level int (0-6)");
    session.cli(cli);

    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) // Indicates a command line error
        return returnCode;

    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));

    return session.run();
}

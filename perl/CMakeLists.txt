cmake_minimum_required(VERSION 2.8)

set(EXECUTABLE_NAME breakdancer-max)
set(SUPPORT_LIBDIR lib/breakdancer-max${EXE_VERSION_SUFFIX})
install(
    FILES AlnParser.pm Poisson.pm bam2cfg.pl
    DESTINATION ${SUPPORT_LIBDIR}
    PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    ) 

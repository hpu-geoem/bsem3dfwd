SET(bsem3dfwd_SOURCES em_ctx.cc em_utils.cc em_io.cc em_la.cc em_response.cc em_fwd.cc em_main.cc)

ADD_EXECUTABLE(bsem3dfwd ${bsem3dfwd_SOURCES})
DEAL_II_SETUP_TARGET(bsem3dfwd)

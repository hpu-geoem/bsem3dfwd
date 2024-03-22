## ----------------------------------------------------------------------------
##
## Copyright (C) 2023 by Ce Qin <ce.qin@hpu.edu.cn>
##
## This file is part of the BSEM3DFWD program.
##
## This code is licensed under the MIT license found in the
## LICENSE file in the root directory of this source tree.
##
## ----------------------------------------------------------------------------

SET(bsem3dfwd_SOURCES em_ctx.cc em_utils.cc em_io.cc em_la.cc em_response.cc em_fwd.cc em_main.cc)

ADD_EXECUTABLE(bsem3dfwd ${bsem3dfwd_SOURCES})
DEAL_II_SETUP_TARGET(bsem3dfwd)

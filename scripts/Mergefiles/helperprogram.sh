#!/bin/env bash
#Connor Depies August, 25, 2017
# This file is for use with makeplinkb.sh only
#list all missing positions
diff --suppress-common-lines -y ffstt.list NAM.list| cut -f 1,2 >missing_both

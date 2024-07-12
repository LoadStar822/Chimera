/*
 * -----------------------------------------------------------------------------
 * Filename:      Chimera.h
 *
 * Author:        Qinzhong Tian
 *
 * Email:         tianqinzhong@qq.com
 *
 * Created Date:  2024-07-09
 *
 * Last Modified: 2024-07-09
 *
 * Description:
 *  This is a simple C++ program that outputs "Hello, World!".
 *
 * Version:
 *  1.0
 * -----------------------------------------------------------------------------
 */
#pragma once

#include <iostream>
#include "kvec.h"
#include "khash.h"
#include "kstring.h"
#include <assert.h>        
#include <math.h>          
#include <iostream>        
#include <vector>          
#include "build/filter/cuckoofilter.h"  
#include <zlib.h>
#include "kseq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include "config.hpp"
#include <string.h>

KSEQ_INIT(gzFile, gzread)

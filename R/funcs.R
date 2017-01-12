# Replacement Bootstrap
#
# Generate a single Bootstrap series according to the sampling-based
# replacement priciple.
#
# Parameters:
# --------------------------------
# seed: Random Seed
# replacement_percentage: The percentage of replacements to perform on the
#                         original series to produce the Bootstrap.
#
# The following algorithm is implemented:
# @inproceedings{sani2015replacement,
# title={The replacement bootstrap for dependent data},
# author={Sani, Amir and Lazaric, Alessandro and Ryabko, Daniil},
# booktitle={2015 IEEE International Symposium on Information Theory (ISIT)},
# pages={1194--1198},
# year={2015},
# organization={IEEE}
# }


.onLoad=function(libname,pkgname){
  lz4Ver = tryCatch( gsub("[^0-9_.-]", "", packageVersion("lz4"), fixed = FALSE),error=function(e) return(-1) );

  if(lz4Ver<0){
    stop("Package lz4 is not installed. You need to install it from github.com as follows:\n\ndevtools::install_github('bwlewis/lz4')\n\nprovided that you have already install devtools.")
  }
}

# Functions
bin_to_float = function(b){
  # Convert binary to float
  return(readBin(packBits(b),"numeric"))
}

float_to_bin = function(f){
  # Convert float to binary
  return(rawToBits(writeBin(f,raw(8))))
}

replacementBootstrap = function(x,
                                seed=NA,
                                replacementFraction=1.0){
  eval(parse(text="library(lz4)"))
  lzCompressLocal=eval(parse(text="lz4::lzCompress"))



  # Check that replacement percentage is less than or equal to 1.0
  stopifnot(replacementFraction<=1)

  # Get original real numbered series length
  if(zoo::is.zoo(x)){
    x=as.numeric(as.vector(zoo::coredata(x)))       #not xts
  }
  x=x[is.finite(x)]   #removes NA and Inf values

  L = length(x)

  # A bit/binary 0-1 representation of the original series
  bit_series = lapply(x,float_to_bin)

  # Construct the conditional distribution
  #
  # appends the same vector, for the sake of compression later in the loop
  #
  # The bit representation concatenated onto the bit representation to provide
  # a conditional distribution for the compressor. All replacements occur on
  # the first half of the conditional distribution list.
  conditional_binary = c(bit_series,bit_series)

  # Set the seed
  if(!is.na(seed)){
    set.seed(seed)
  }

  # The replacement index
  # A permutation of size L
  replacement_points = sample(seq_along(x),size=L)

  # Preallocate containers
  probabilities=rep(NA,L)
  bootstrap_series=rep(NA,L)

  # Replacement Bootstrap Loop
  # Compute the probability distribution for a specific replacement index over the alphabet of possible replacements realized in the given series
  for(replacement_idx in replacement_points){
    for(i in seq_along(bit_series)){
      # 1. Replace the replacement_idx point in the series with each of the
      # bit representations of the original series.
      conditional_binary[[replacement_idx]] = bit_series[[i]]

      # 2. Compress each new replaced series and invert to get
      # approximate probability of that bit representation falling in
      # replacement_idx according to the process that generated the series.
      probabilities[i] = 1.0/length(lzCompressLocal(unlist(conditional_binary)))
    }
    # Normalize the probabilities
    probabilities = probabilities/sum(probabilities)

    # Draw a single replacement according to the replacement distribution.
    replacement_choice = sample(L, 1, prob = probabilities)

    # Replace the existing value in the original series with the replacement
    # value drawn from the replacement distribution.
    bootstrap_series[replacement_idx] = x[replacement_choice]

    # Replace the bit representation in replacement_idx with the bit
    # representation of the replacement value.
    conditional_binary[replacement_idx] = bit_series[replacement_choice]
  }
  return(bootstrap_series)
}


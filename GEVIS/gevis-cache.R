# Set the cache directory and maximum cache age in seconds
cache_dir <- "/tmp/ocpu-temp/.cache/R/R.cache"
max_age <- 2000  # 1 hour
check_interval <- 600  # Time between checks in seconds (e.g., 10 minutes)

# Helper function to get the current time in seconds
current_time <- function() {
  as.numeric(Sys.time())
}

# Function to clean up expired cache files
cleanup_expired_cache <- function() {
  # Get all files in the cache directory, ignoring README.txt
  all_cache_files <- list.files(cache_dir, recursive = TRUE, full.names = TRUE)
  cache_files <- all_cache_files[!grepl("README.txt$", all_cache_files)]

  # Print all cache files before cleanup
  message("Listing all cache files:")
  print(cache_files)

  for (file in cache_files) {
    # Get file modification time and convert it to seconds
    file_info <- file.info(file)
    cache_age <- current_time() - as.numeric(file_info$mtime)

    # Check if file is older than max_age
    if (cache_age > max_age) {
      message("Removing expired cache file: ", file)
      file.remove(file)
    }
  }
}

# Continuously run the cache cleanup at specified intervals
while (TRUE) {
  message("Running cache cleanup...")
  cleanup_expired_cache()
  Sys.sleep(check_interval)  # Pause between cleanups
}

# Simple try operator - returns default value if expression fails
`%try%` <- function(expr, fallback_expr = NULL) {
      result <- NULL
      tryCatch({
            result <- eval(substitute(expr), envir = parent.frame())
            result
      }, error = function(e) {
            fallback_expr
      })
}

# Retry operator - attempts multiple times before using fallback
`%retry%` <- function(expr, fallback_expr = NULL, attempts = 3, delay = 1) {

      for (i in 1:attempts) {
            result <- tryCatch({
                  eval(substitute(expr), envir = parent.frame())
            }, error = function(e) {
                  message(sprintf("Attempt %d failed: %s", i, e$message))
                  if (i == attempts) return(NULL)  # Flag for final attempt failed
                  Sys.sleep(delay)
                  return(NULL)  
            })
            
            if (!is.null(result)) return(result)  # Return successful result
      }
      
      return(fallback_expr)  # Return fallback value if all attempts fail
}

test_constructor <- function(){
  cfg <- cobindRConfiguration() # just check this does not throw an error
  checkException(cfg <- cobindRConfiguration('nonexistingFile.txt'))
}

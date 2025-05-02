isWindows <- function()
{
    ret <- FALSE
    sysinf <- Sys.info()
    if (sysinf["sysname"] == 'Windows')
    {
        ret <- TRUE
    }
    
    return(ret)
}

getDllName <- function(dll_base_name)
{
    dll_name =  paste(dll_base_name,".so",sep="")
                      
    if(isWindows())
    {
        dll_name =  paste(dll_base_name,".dll",sep="")
    } 
    return(dll_name)
}

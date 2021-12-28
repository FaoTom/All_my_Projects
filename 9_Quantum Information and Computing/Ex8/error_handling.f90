MODULE error_handling

    ! MODULE error_handling
    ! -> This module error_handling can be used to debug the code
    !    and deal with errors of different natures that may arise
    !    in the execution of the code.
    ! 
    ! -> To activate debug mode add '-debug' flag when launching
    !    the executable and call checkpoint SUBROUTINE when needed.
    !
    ! =================================================================
    !
    ! FUNCTIONS AND SUBROUTINES:
    ! - SUBROUTINE checkpoint(debug, msg, input)
    !

CONTAINS

SUBROUTINE checkpoint(debug, msg, type, input)

    ! SUBROUTINE checkpoint
    !  -> This function prints a message (if passed) of different
    !     nature and a variable (if passed) with its correspondent TYPE.
    !
    ! -> To activate debug mode add '-debug' flag when launching
    !    the executable and call checkpoint SUBROUTINE when needed.
    !
    ! =================================================================
    !
    ! PARAMETERS:
    ! - debug: boolean value for entering debug mode (passed via command line)
    ! - msg:   if debug, this optional message is printed
    ! - type:  if debug, this optional arguments can be 'e' for ERROR, and 'd' for DEBUG
    ! - input: if debug, this optional input of unspecified TYPE is printed
    !

IMPLICIT NONE

LOGICAL,      INTENT(IN)           :: debug
!optional arguments
CHARACTER(*), INTENT(IN), OPTIONAL :: msg, type
CLASS(*),     INTENT(IN), OPTIONAL :: input

!if debug variable is TRUE
IF(debug.EQV..TRUE.) THEN
    !at first, check if a message is passed
    IF(PRESENT(msg)) THEN
        !then, look whether the type of the message is specified
        IF (PRESENT(type)) THEN
            IF(type.EQ.'e') THEN
                !for ERROR, write in bold red [ERROR] before the message
                PRINT *, achar(27)//"[1;31m"//"[ERROR]"//achar(27)//"[0m "//msg
            ELSE IF(type.EQ.'d') THEN 
                !for ERROR, write in bold red [ERROR] before the message
                PRINT *, achar(27)//"[1;32m"//"[DEBUG]"//achar(27)//"[0m "//msg
            !for a different type print the message alone
            ELSE
                PRINT *, msg
        ENDIF
        !if type is not specified print the message alone
        ELSE
            PRINT *, msg
        ENDIF
    ENDIF

    !if an input variable is passed
    IF(PRESENT(input)) then
        !print the variable and its TYPE accordingly
        SELECT TYPE (input)
        TYPE IS (LOGICAL)
        WRITE (*, '(L)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> LOGICAL variable"
        TYPE IS (INTEGER(2))
        WRITE (*, '(I2)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> INTEGER(2) variable"
        TYPE IS (INTEGER(4))
        WRITE (*, '(I4)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> INTEGER(4) variable"
        TYPE IS (INTEGER(8))
        WRITE (*, '(I8)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> INTEGER(8) variable"
        TYPE IS (REAL(4))
        WRITE (*, '(I8)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> REAL(4) variable"
        TYPE IS (REAL(8))
        WRITE (*, '(I8)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> REAL(8) variable"
        TYPE IS (COMPLEX(4))
        WRITE (*, '(I8)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> COMPLEX(4) variable"
        TYPE IS (COMPLEX(8))
        WRITE (*, '(I8)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> COMPLEX(8) variable"
        TYPE IS (COMPLEX(16))
        WRITE (*, '(I8)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> COMPLEX(16) variable"
        TYPE IS (CHARACTER(*))
        WRITE (*, '(I8)', ADVANCE = 'no') input
        WRITE (*, '(A)') " -> CHARACTER(*) variable"
    END SELECT
    ENDIF 
ENDIF

END SUBROUTINE checkpoint

END MODULE error_handling
# Run a command detached from the terminal session (survives the session
# ending): python3 detach.py <logfile> <command> [args...]
import os, sys
log = sys.argv[1]
if os.fork() == 0:
    os.setsid()
    fd = os.open(log, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
    os.dup2(fd, 1); os.dup2(fd, 2)
    os.environ.update(VECLIB_MAXIMUM_THREADS="1", OMP_NUM_THREADS="1",
                      OPENBLAS_NUM_THREADS="1")
    os.execvp(sys.argv[2], sys.argv[2:])

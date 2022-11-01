#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

void check_err_open_lseek(int result) {
  if(result == -1) {
    perror("Error opening file/lseek for reading");
    exit(EXIT_FAILURE);
  }
}

void check_fail_mmap(double* map, int fd) {
  if(map == MAP_FAILED) {
    close(fd);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
  }
}

void check_err_write(int result, int fd) {
  if(result != 1) {
    close(fd);
    perror("Error writing last byte of the file");
    exit(EXIT_FAILURE);
  }
}


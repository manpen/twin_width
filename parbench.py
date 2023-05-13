import os
import subprocess
import signal
import sys
import time

def run_binary_with_files(dir1_path, binary_path, dir2_path, timeout, max_processes):
    # Validate the paths
    dir1_path = os.path.join(os.getcwd(), dir1_path)
    binary_path = os.path.join(os.getcwd(), binary_path)
    dir2_path = os.path.join(os.getcwd(), dir2_path)

    if not os.path.isdir(dir1_path):
        print(f"Error: '{dir1_path}' is not a valid directory path.")
        return
    if not os.path.isfile(binary_path):
        print(f"Error: '{binary_path}' is not a valid file path.")
        return
    if not os.path.isdir(dir2_path):
        print(f"Error: '{dir2_path}' is not a valid directory path.")
        return

    # Initialize counters and process list
    total_files = 0
    processed_files = 0
    processes = []

    # Iterate over files in dir1
    for filename in os.listdir(dir1_path):
        file_path = os.path.join(dir1_path, filename)
        if os.path.isfile(file_path):
            total_files += 1

    # Iterate over files in dir1
    for filename in os.listdir(dir1_path):
        file_path = os.path.join(dir1_path, filename)
        if os.path.isfile(file_path):
            # Create the output file paths in dir2
            stdout_path = os.path.join(dir2_path, f"{filename}.out")
            stderr_path = os.path.join(dir2_path, f"{filename}.err")

            # Run the binary as a background process
            cmd = [binary_path]
            with open(file_path, "r") as input_file, \
                    open(stdout_path, "w") as stdout_file, \
                    open(stderr_path, "w") as stderr_file:
                process = subprocess.Popen(cmd, stdin=input_file, stdout=stdout_file, stderr=stderr_file, preexec_fn=os.setsid)
                processes.append(process)

                # Increment the processed files counter
                processed_files += 1

                # Check if the maximum processes limit is reached or if all files are processed
                if processed_files % max_processes == 0 or processed_files == total_files:
                    # Wait for the specified timeout
                    time.sleep(timeout)

                    # Terminate each process if it's still running
                    for process in processes:
                        process.terminate()

                    # Wait for 5 seconds
                    time.sleep(5)

                    # Kill each process if it's still running
                    for process in processes:
                        process.kill()

                    # Print process information
                    print(f"Termination sequence initiated for batch {processed_files // max_processes + 1}")

                    # Clear the processes list
                    processes = []

    # Print overall process information
    print(f"All processes terminated.")

if __name__ == "__main__":
    # Check if all five arguments are provided
    if len(sys.argv) != 6:
        print("Usage: python parbench.py <binary_path> <input_dir> <output_dir> <timeout> <max_processes>")
    else:
        # Extract the arguments
        binary_path = sys.argv[1]
        dir1_path = sys.argv[2]
        dir2_path = sys.argv[3]
        timeout = int(sys.argv[4])
        max_processes = int(sys.argv[5])

        # Run the binary with files
        run_binary_with_files(dir1_path, binary_path, dir2_path, timeout, max_processes)

import os
import subprocess
import signal
import sys
import time
import csv
import argparse
import pace_verifier

def validate_directories(dir1_path, binary_path, dir2_path):
    # Validate the paths
    dir1_path = os.path.join(os.getcwd(), dir1_path)
    binary_path = os.path.join(os.getcwd(), binary_path)
    dir2_path = os.path.join(os.getcwd(), dir2_path)

    if not os.path.isdir(dir1_path):
        print(f"Error: '{dir1_path}' is not a valid directory path.")
        return False
    if not os.path.isfile(binary_path):
        print(f"Error: '{binary_path}' is not a valid file path.")
        return False
    if not os.path.isdir(dir2_path):
        print(f"Error: '{dir2_path}' is not a valid directory path.")
        return False

    return True

def run_binary_with_files(input_dir, binary_path, output_dir, timeout, num_threads, kill_buffer):
    # Initialize counters and process list
    total_files = 0
    processed_files = 0
    processes = []

    # Create a stack of all files
    file_stack = []
    for filename in os.listdir(input_dir):
        file_path = os.path.join(input_dir, filename)
        if os.path.isfile(file_path):
            file_stack.append(file_path)
            total_files += 1

    with open(os.path.join(output_dir, "execution_stats.csv"), "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["name", "time"])

        # Loop until all files are processed and no active processes remain
        while file_stack or processes:
            # Check if any processes have terminated
            for i, (process, start_time, sigterm_flag) in enumerate(processes):
                if process.poll() is not None:
                    end_time = time.time()
                    runtime = end_time - start_time
                    writer.writerow([os.path.basename(file_path), runtime])
                    processes.pop(i)
                    processed_files += 1
                    print(f"Process terminated with runtime: {runtime:.2f} seconds")

            # Replace terminated processes with new ones (if possible)
            while len(processes) < num_threads and file_stack:
                file_path = file_stack.pop()
                stdout_path = os.path.join(output_dir, f"{os.path.basename(file_path)}.out")
                stderr_path = os.path.join(output_dir, f"{os.path.basename(file_path)}.err")

                cmd = [binary_path]
                with open(file_path, "r") as input_file, \
                        open(stdout_path, "w") as stdout_file, \
                        open(stderr_path, "w") as stderr_file:
                    process_start_time = time.time()
                    process = subprocess.Popen(cmd, stdin=input_file, stdout=stdout_file, stderr=stderr_file,
                                               preexec_fn=os.setsid)
                    processes.append((process, process_start_time, False))

            # Check timeout and terminate or kill processes
            for i in range(len(processes)):
                process, start_time, sigterm_flag = processes[i]
                elapsed_time = time.time() - start_time
                if elapsed_time > timeout:
                    if sigterm_flag and elapsed_time > timeout + kill_buffer:
                        process.kill()
                    elif not sigterm_flag:
                        process.terminate()
                        processes[i] = (process, start_time, True)

            # Sleep for a second before checking process status again
            time.sleep(1)

    # Print overall process information
    print(f"All processes terminated.")

# Process files in dir1_path and dir2_path and create a results.csv file
def process_files(dir1_path, dir2_path):
    results = []
    for filename in os.listdir(dir1_path):
        file_path1 = os.path.join(dir1_path, filename)
        file_path2 = os.path.join(dir2_path, f"{filename}.out")
        if os.path.isfile(file_path1) and os.path.isfile(file_path2):
            try:
                graph = pace_verifier.read_graph(file_path1)
                sequence = pace_verifier.read_sequence(file_path2)
                result = pace_verifier.check_sequence(graph, sequence)
                results.append((filename, result))
            except Exception as e:
                results.append((filename, "EXCEPTION"))
        else:
            results.append((filename, "MISSING"))

    # Write results to CSV file
    csv_path = os.path.join(dir2_path, "results.csv")
    with open(csv_path, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["name", "tww"])
        writer.writerows(results)

# Main function
def main(args):
    input_dir = args.input
    binary_path = args.binary
    output_dir = args.output
    timeout = args.timeout
    num_threads = args.num_threads
    verify_only = args.verify
    kill_buffer = args.kill_buffer

    # Skip the execution of the solver if only the verifier should be called
    if not verify_only:
        run_binary_with_files(input_dir, binary_path, output_dir, timeout, num_threads, kill_buffer)

    process_files(input_dir, output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input directory")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("-t", "--timeout", type=int, default=300, help="Timeout in seconds (default: 300)")
    parser.add_argument("-n", "--num_threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("-k", "--kill_buffer", type=int, default=5, help="Kill buffer in seconds (default: 5)")
    parser.add_argument("--verify", action="store_true", default=False, help="Verify only, no binary execution")
    parser.add_argument("--binary", type=str, default=False, required=True, help="Binary path of solver")
    args = parser.parse_args()
    main(args)

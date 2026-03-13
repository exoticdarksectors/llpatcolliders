import os
import shlex
import subprocess
from pathlib import Path


MG5_IMAGE = "mg5-hnl"


def docker_repo_path(project_root, host_path):
    """Map a repo-local host path into the container worktree."""
    return Path("/work") / Path(host_path).relative_to(project_root)


def docker_command(project_root, command, cwd=None):
    """Build a Docker command that executes inside the repo checkout."""
    if cwd is None:
        cwd = project_root

    cmd = ["docker", "run", "--rm"]
    if hasattr(os, "getuid") and hasattr(os, "getgid"):
        cmd += ["--user", f"{os.getuid()}:{os.getgid()}"]
    cmd += [
        "-e", "HOME=/tmp",
        "-v", f"{project_root}:/work",
        "-w", str(docker_repo_path(project_root, cwd)),
        MG5_IMAGE,
        "/bin/bash", "-lc", command,
    ]
    return cmd


def run_mg5_command_file(project_root, cmd_file, log_file, timeout):
    """Run an MG5 command file inside the Docker image."""
    cmd = docker_command(
        project_root,
        f"set -euo pipefail; mg5_aMC {shlex.quote(str(docker_repo_path(project_root, cmd_file)))}",
    )
    with open(log_file, "w") as log:
        return subprocess.run(
            cmd, stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
        )


def run_generate_events(project_root, work_subdir, log_file, timeout, nb_core=1):
    """Run the generated MG5 process inside the Docker image."""
    command = [
        "python3", "bin/generate_events",
        "-f", "--laststep=parton",
    ]
    if nb_core > 1:
        command += ["--multicore", f"--nb_core={nb_core}"]
    else:
        command += ["--nb_core=1"]

    shell_command = "set -euo pipefail; " + " ".join(shlex.quote(part) for part in command)
    cmd = docker_command(project_root, shell_command, cwd=work_subdir)
    with open(log_file, "w") as log:
        return subprocess.run(
            cmd, stdout=log, stderr=subprocess.STDOUT, timeout=timeout,
        )

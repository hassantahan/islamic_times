# Release Checklist

1. Review the changelog and the upgrade implications of changed defaults.
2. Run CI on the release pull request. Manually run **Build and Publish Wheels**
   on the release branch to check the distributions before publishing.
3. Verify CPython 3.10 through 3.13 wheels for Windows (x86 and x86-64), Linux
   (x86-64), and macOS (Intel and Apple Silicon). Check the wheel and source
   distribution metadata. Install the candidate wheel and source distribution
   in clean environments and check that both import `islamic_times.astro_core`.
   Verify the `[map]`
   extra and the `islamic-times-map` command.
4. After merge, confirm CI passes on the target commit on `master`.
5. Tag that commit `vX.Y.Z` and push the tag. The tag sets the distribution
   version and triggers the PyPI publish job.
6. Verify the new version and files on PyPI, then install the exact version in
   a clean environment. Verify dependent applications against that install
   before deploying them.

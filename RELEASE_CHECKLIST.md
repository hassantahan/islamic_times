# Release Checklist

1. Ensure CI is green on `master` for the target commit.
2. Confirm versioning/tag plan (semantic version and release notes).
3. Verify core install and optional mapping install paths:
   - `pip install islamic_times`
   - `pip install "islamic_times[map]"`
4. Confirm wheel and sdist builds succeed in GitHub Actions.
5. Create and push tag `vX.Y.Z` to trigger publish workflow.
6. Verify package artifacts and metadata on PyPI.

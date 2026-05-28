#!/usr/bin/env bash
# push-both.sh — push the current branch to both git remotes.
#
# This repo has two remotes:
#   - origin     = https://github.com/genedance/speed-dada.git    (public)
#   - wurstator  = git@…:wurstator/speed-dada.git                 (private)
#
# Both share an identical `main` branch and identical tags. Claude-only
# state lives on the orphan `claude-state` branch in `wurstator` and is
# never pushed to `origin` (see `wurstator/claude-state` → README.md for
# the rationale). To refresh that branch, run `tools/sync-claude-state.sh`
# from the `claude-state` worktree, not this script.
#
# Usage:
#   bash tools/push-both.sh              # pushes current branch + matching tags
#   bash tools/push-both.sh --tags-only  # pushes only tags

set -euo pipefail

BRANCH="$(git branch --show-current)"
if [ -z "$BRANCH" ]; then
    echo "ERROR: detached HEAD — check out a branch first" >&2
    exit 1
fi
if [ "$BRANCH" = "claude-state" ]; then
    echo "ERROR: refusing to push claude-state to both remotes." >&2
    echo "       claude-state lives only in wurstator. Push it with:" >&2
    echo "         git push wurstator claude-state" >&2
    exit 1
fi

case "${1:-}" in
    --tags-only)
        git push origin --tags
        git push wurstator --tags
        ;;
    "")
        git push origin "$BRANCH"
        git push wurstator "$BRANCH"
        git push origin --tags
        git push wurstator --tags
        ;;
    *)
        echo "Usage: $0 [--tags-only]" >&2
        exit 1
        ;;
esac

echo "push-both: $BRANCH + tags pushed to origin (genedance) and wurstator."

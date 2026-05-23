#!/usr/bin/env bash
# Fires before every Bash command. Blocks destructive patterns; logs all commands.
set -euo pipefail

CMD="${CLAUDE_TOOL_INPUT:-}"
LOG="/home/fromage/dada2_rust/.claude/hooks/command.log"

echo "$(date -u +%FT%TZ) PRE  $CMD" >> "$LOG"

# Hard-block patterns
if echo "$CMD" | grep -qE 'rm\s+-rf\s+/|git\s+push\s+--force|curl\s|wget\s'; then
    echo "BLOCKED: dangerous command pattern detected" >&2
    exit 2
fi

exit 0

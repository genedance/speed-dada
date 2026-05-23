#!/usr/bin/env bash
# Fires after every file write. Auto-formats Rust files; warns if > 500 lines.
set -euo pipefail

FILE="${CLAUDE_TOOL_RESULT_FILE:-}"
LOG="/home/fromage/dada2_rust/.claude/hooks/writes.log"

echo "$(date -u +%FT%TZ) WRITE $FILE" >> "$LOG"

if [[ "$FILE" == *.rs ]]; then
    rustfmt "$FILE" 2>/dev/null || true
    LINES=$(wc -l < "$FILE")
    if (( LINES > 500 )); then
        echo "WARNING: $FILE has $LINES lines (limit 500). Consider splitting." >&2
    fi
fi

exit 0

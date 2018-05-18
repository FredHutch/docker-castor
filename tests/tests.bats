#!/usr/bin/env bats

@test "HSP wrapper script" {
  v="$(castor_hidden_state_prediction.R -h 2>&1 || true )"
  [[ "$v" =~ "Perform hidden state prediction" ]]
}

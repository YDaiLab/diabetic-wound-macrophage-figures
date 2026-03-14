# --- monkey patch ggh4x internal validator (SESSION ONLY) ---
ggh4x_ns <- asNamespace("ggh4x")

# save original
old_validate <- get("validate_element_list", envir = ggh4x_ns)

# replace with a permissive version (just returns input unchanged)
unlockBinding("validate_element_list", ggh4x_ns)
assign(
  "validate_element_list",
  function(x, ...) x,
  envir = ggh4x_ns
)
lockBinding("validate_element_list", ggh4x_ns)
# ------------------------------------------------------------

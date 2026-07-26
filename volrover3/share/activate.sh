# VolRover3 environment — source this:   source <prefix>/bin/activate
#
# Puts the installed `volrover3` (and any other bundled commands) on PATH and sets
# the library paths its shared libs + embedded Python need, so they run with NO
# further configuration. The interpreter self-locates its CPython home and the
# vrhost module from the binary's own location, so this script only wires up the
# dynamic-library search path. Run `deactivate` to undo.
#
# Works in bash and zsh.

# --- resolve this bundle's prefix from the script's own location (<prefix>/bin/activate) ---
if [ -n "${BASH_SOURCE:-}" ]; then
  _vr_self="${BASH_SOURCE[0]}"
elif [ -n "${ZSH_VERSION:-}" ]; then
  _vr_self="${(%):-%x}"
else
  _vr_self="$0"
fi
_vr_bin="$(cd "$(dirname "$_vr_self")" >/dev/null 2>&1 && pwd)"
VOLROVER3_PREFIX="$(dirname "$_vr_bin")"

# --- remember what we override, so `deactivate` can restore it exactly ---
_VR_OLD_PATH="$PATH"
_VR_OLD_LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"
_VR_OLD_QT_PLUGIN_PATH="${QT_PLUGIN_PATH:-}"
_VR_OLD_PS1="${PS1:-}"
_VR_HAD_LDP="${LD_LIBRARY_PATH+set}"
_VR_HAD_QPP="${QT_PLUGIN_PATH+set}"

export PATH="$_vr_bin:$PATH"

# The bundle's own libs (cvc, cvcGL, VTK, Qt, ...) plus numpy's vendored BLAS.
_vr_ldp="$VOLROVER3_PREFIX/lib"
for _np in "$VOLROVER3_PREFIX"/lib/python*/site-packages/numpy.libs; do
  [ -d "$_np" ] && _vr_ldp="$_vr_ldp:$_np"
done
export LD_LIBRARY_PATH="$_vr_ldp${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# Qt plugins, if this bundle ships them (else Qt finds them next to its libs).
if [ -d "$VOLROVER3_PREFIX/plugins" ]; then
  export QT_PLUGIN_PATH="$VOLROVER3_PREFIX/plugins${QT_PLUGIN_PATH:+:$QT_PLUGIN_PATH}"
fi

export VOLROVER3_PREFIX

deactivate() {
  export PATH="$_VR_OLD_PATH"
  if [ "$_VR_HAD_LDP" = set ]; then export LD_LIBRARY_PATH="$_VR_OLD_LD_LIBRARY_PATH"; else unset LD_LIBRARY_PATH; fi
  if [ "$_VR_HAD_QPP" = set ]; then export QT_PLUGIN_PATH="$_VR_OLD_QT_PLUGIN_PATH"; else unset QT_PLUGIN_PATH; fi
  [ -n "$_VR_OLD_PS1" ] && PS1="$_VR_OLD_PS1"
  unset VOLROVER3_PREFIX _vr_self _vr_bin _vr_ldp _np
  unset _VR_OLD_PATH _VR_OLD_LD_LIBRARY_PATH _VR_OLD_QT_PLUGIN_PATH _VR_OLD_PS1 _VR_HAD_LDP _VR_HAD_QPP
  unset -f deactivate
}

case "${PS1:-}" in
  *VolRover3*) ;;
  *) PS1="(VolRover3) ${PS1:-}" ;;
esac

echo "VolRover3 environment activated: $VOLROVER3_PREFIX"
echo "  installed commands are on PATH (e.g. 'volrover3'); example scripts in share/volrover3/examples; run 'deactivate' to undo."

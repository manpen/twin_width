pub fn checkIfIsCompatibleInteger(comptime T: type) bool {
    comptime if (T != u8 and T != u16 and T != u32) {
        return false;
    };
    return true;
}

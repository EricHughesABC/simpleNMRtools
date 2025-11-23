function testDir() {
    var newDir = Dir.appDataCommon();
    print("Dir.appDataCommon():\n\t" + newDir);
    print("Dir.appDataUser ( ):\n\t" + Dir.appDataUser());
    print("Dir.application ( ):\n\t" + Dir.application());
    print("Dir.searchPaths ( String aPrefix ):\n\t" + Dir.searchPaths(".qs"));
}
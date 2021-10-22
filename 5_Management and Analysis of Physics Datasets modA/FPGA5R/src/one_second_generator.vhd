LIBRARY ieee;
USE ieee.std_logic_1164.ALL;
USE ieee.numeric_std.ALL;

----NOT USED IN PROJECT
ENTITY one_second_generator IS

  PORT (
    clock : IN STD_LOGIC;
    one_second_out : OUT STD_LOGIC);

END ENTITY one_second_generator;
ARCHITECTURE rtl OF one_second_generator IS
  SIGNAL counter : unsigned(31 DOWNTO 0) := (OTHERS => '0');
  CONSTANT divisor : unsigned(31 DOWNTO 0) := to_unsigned(100000000, 32);
BEGIN -- architecture rtl
  main : PROCESS (clock) IS
  BEGIN -- process main
    IF rising_edge(clock) THEN -- rising clock edge
      counter <= counter + 1;
      IF counter = divisor THEN
        one_second_out <= '1';
        counter <= (OTHERS => '0');
      ELSE
        one_second_out <= '0';
      END IF;
    END IF;
  END PROCESS main;

END ARCHITECTURE rtl;
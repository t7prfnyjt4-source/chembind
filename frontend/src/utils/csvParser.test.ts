import { describe, it, expect } from "vitest";
import { parseSmilesCsv } from "./csvParser";

describe("parseSmilesCsv - sanity + edge cases", () => {
  it("parses with header", () => {
    const csv = `smiles
CCO
CCC`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
    expect(rows[0]).toEqual({ rowId: "row_2", smiles: "CCO" });
    expect(rows[1]).toEqual({ rowId: "row_3", smiles: "CCC" });
  });

  it("parses without header", () => {
    const csv = `CCO
CCC`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
    expect(rows[0].smiles).toBe("CCO");
    expect(rows[1].smiles).toBe("CCC");
  });

  it("ignores empty lines", () => {
    const csv = `smiles

CCO

CCC

`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
  });

  it("handles uppercase header", () => {
    const csv = `SMILES
CCO
CCC`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
  });

  it("handles header with additional columns", () => {
    const csv = `id,smiles,name
1,CCO,ethanol
2,CCC,propane`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
    expect(rows[0].smiles).toBe("CCO");
    expect(rows[1].smiles).toBe("CCC");
  });

  it("handles quoted values with commas", () => {
    const csv = `smiles
"CCO,CO"
"CCC"`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
    expect(rows[0].smiles).toBe("CCO,CO");
  });

  it("filters rows where smiles column is empty", () => {
    const csv = `smiles
CCO

CCC
`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
  });

  it("returns empty array for empty file", () => {
    expect(parseSmilesCsv("")).toEqual([]);
    expect(parseSmilesCsv("   \n \n")).toEqual([]);
  });

  it("handles BOM at file start", () => {
    const csv = `\uFEFFsmiles
CCO
CCC`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(2);
  });

  it("preserves special SMILES characters", () => {
    const csv = `smiles
C1=CC=CC=C1
CC(=O)O
[NH4+]`;
    const rows = parseSmilesCsv(csv);

    expect(rows.length).toBe(3);
    expect(rows[2].smiles).toBe("[NH4+]");
  });
});

import { useEffect, useState } from "react";

import {
  type User,
  onAuthStateChanged,
  GoogleAuthProvider,
  signInWithPopup,
  signOut,
} from "firebase/auth";

import { auth } from "./firebase/config";


type Page = "analyze" | "batch" | "history" | "search";

export default function App() {
  const [user, setUser] = useState<User | null>(null);
  const [authLoading, setAuthLoading] = useState(true);
  const [page, setPage] = useState<Page>("analyze");

  useEffect(() => {
    const unsub = onAuthStateChanged(auth, (u) => {
      setUser(u);
      setAuthLoading(false);
    });
    return () => unsub();
  }, []);

  const login = async () => {
    const provider = new GoogleAuthProvider();
    await signInWithPopup(auth, provider);
  };

  const logout = async () => {
    await signOut(auth);
  };

  if (authLoading) return null;

  return (
    <div style={{ padding: 24, fontFamily: "system-ui" }}>
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
          gap: 12,
        }}
      >
        <h1 style={{ margin: 0 }}>Chembind</h1>

        <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
          {user ? (
            <>
              <span style={{ opacity: 0.8, fontSize: 14 }}>
                {user.email ?? user.uid}
              </span>
              <button onClick={logout}>Logout</button>
            </>
          ) : (
            <button onClick={login}>Login</button>
          )}
        </div>
      </div>

      <div style={{ display: "flex", gap: 8, marginTop: 16 }}>
        <button onClick={() => setPage("analyze")}>Analyze</button>
        <button onClick={() => setPage("batch")}>Batch</button>
        <button onClick={() => setPage("history")}>History</button>
        <button onClick={() => setPage("search")}>Search</button>
      </div>

      <div
        style={{
          marginTop: 16,
          padding: 12,
          border: "1px solid #333",
          borderRadius: 8,
        }}
      >
        {page === "analyze" && <div>Analyze page placeholder</div>}
        {page === "batch" && <div>Batch page placeholder</div>}
        {page === "history" && <div>History page placeholder</div>}
        {page === "search" && <div>Search page placeholder</div>}
      </div>
    </div>
  );
}

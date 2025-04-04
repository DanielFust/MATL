//
//  testApp.swift
//  Shared
//
//  Created by Daniel Fust on 3/13/22.
//

import SwiftUI

@main
struct testApp: App {
    var body: some Scene {
        DocumentGroup(newDocument: testDocument()) { file in
            ContentView(document: file.$document)
        }
    }
}
